!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module equation_of_state
  !!< This module contains functions used to evaluate the equation of state.
  use fldebug
  use fields
  use state_module
  use global_parameters, only: OPTION_PATH_LEN
  use spud
  use sediment, only: get_sediment_name, get_nSediments
  use boundary_conditions, only: set_dirichlet_consistent, get_boundary_condition_count
  use petsc_solve_state_module
  use sparsity_patterns_meshes
  use transform_elements, only: transform_to_physical
  use diagnostic_fields, only : calculate_galerkin_projection
  
  implicit none
  
  interface compressible_eos
     module procedure compressible_eos_1mat,&
          compressible_eos_mmat
  end interface compressible_eos

  private
  public :: calculate_perturbation_density, mcD_J_W_F2002, &
            compressible_eos, compressible_material_eos, &
            set_EOS_pressure_and_temperature

contains

  

  subroutine calculate_perturbation_density(state, density, reference_density)
    !!< Calculates the perturbation density (i.e. the reference density is already subtracted)
    !!< of a state with equation_of_state fluids/linear or 
    !!< fluids/ocean_pade_approximation.
    type(state_type), intent(in):: state
    type(scalar_field), intent(inout) :: density
    real, intent(out), optional :: reference_density
    
    type(vector_field), pointer:: u
    type(scalar_field), pointer:: T, S, oldT, oldS, topdis
    type(scalar_field) DeltaT, DeltaS, remapT, remapS, fluidconcentration,&
         & sedimentdensity
    character(len=OPTION_PATH_LEN) option_path, dep_option_path, class_name
    logical, dimension(:), allocatable:: done
    logical include_depth_below
    real T0, S0, gamma, rho_0, salt, temp, dist, dens, theta
    integer, dimension(:), pointer:: density_nodes
    integer ele, i, node, nSediments
    
    ewrite(1,*) 'In calculate_perturbation_density'
    
    u => extract_vector_field(state, "Velocity")
    
    call zero(density)
    
    option_path='/material_phase::'//trim(state%name)//'/equation_of_state/fluids'
    
    call get_option(trim(u%option_path)//'/prognostic/temporal_discretisation/relaxation', &
                    theta, default = 1.0)
    
    rho_0 = 0.0
    
    if (have_option(trim(option_path)//'/linear')) then
    
       option_path=trim(option_path)//'/linear'
       
       if (have_option(trim(option_path)//'/temperature_dependency')) then
          dep_option_path=trim(option_path)//'/temperature_dependency'
          call get_option(trim(dep_option_path)//'/reference_temperature', T0)
          call get_option(trim(dep_option_path)//'/thermal_expansion_coefficient', gamma)
          T => extract_scalar_field(state, "Temperature")
          oldT => extract_scalar_field(state, "OldTemperature")
          call allocate(deltaT, density%mesh, "DeltaT")
          call allocate(remapT, density%mesh, "RemapT")
          
          ! deltaT=theta*T+(1-theta)*oldT-T0
          call remap_field(T, remapT)
          call set(deltaT, remapT)
          call scale(deltaT, theta)
          
          call remap_field(oldT, remapT)
          call addto(deltaT, remapT, 1.0-theta)
          call addto(deltaT, -T0)
          ! density=density-gamma*deltaT
          call addto(density, deltaT, scale=-gamma)
          call deallocate(deltaT)
          call deallocate(remapT)
       end if
       
       if (have_option(trim(option_path)//'/salinity_dependency')) then
          dep_option_path=trim(option_path)//'/salinity_dependency'
          call get_option(trim(dep_option_path)//'/reference_salinity', S0)
          call get_option(trim(dep_option_path)//'/saline_contraction_coefficient', gamma)
          S => extract_scalar_field(state, "Salinity")
          oldS => extract_scalar_field(state, "OldSalinity")
          call allocate(deltaS, density%mesh, "DeltaS")
          call allocate(remapS, density%mesh, "RemapS")
          
          ! deltaS=theta*S+(1-theta)*oldS-S0
          call remap_field(S, remapS)
          call set(deltaS, remapS)
          call scale(deltaS, theta)
          
          call remap_field(oldS, remapS)
          call addto(deltaS, remapS, 1.0-theta)
          call addto(deltaS, -S0)
          ! density=density+gamma*deltaS
          call addto(density, deltaS, scale=gamma)
          call deallocate(deltaS)
          call deallocate(remapS)
       end if
       
       call get_option(trim(option_path)//'/reference_density', rho_0)
       call scale(density, rho_0)
       
    elseif (have_option(trim(option_path)//'/ocean_pade_approximation')) then
      
       option_path=trim(option_path)//'/ocean_pade_approximation'
       
       include_depth_below=have_option(trim(option_path)//'/include_depth_below_surface')
       
       T => extract_scalar_field(state, "Temperature")
       oldT => extract_scalar_field(state, "OldTemperature")
       S => extract_scalar_field(state, "Salinity")
       oldS => extract_scalar_field(state, "OldSalinity")
       if (include_depth_below) then
          topdis => extract_scalar_field(state, "DistanceToTop")
       endif
      
       allocate( done(1:node_count(density)) )
       done=.false.
       
       do ele=1, element_count(density)

          density_nodes => ele_nodes(density, ele)

          do i=1,size(density_nodes)
             node=density_nodes(i)
             ! In the continuous case ensure we only do each calculation once.
             if (done(node)) cycle
             done(node)=.true.
            
             salt=theta*node_val(S, node)+(1-theta)*node_val(oldS, node)
             temp=node_val(T, node)+(1-theta)*node_val(oldT, node)
             if (include_depth_below) then
                dist=node_val(topdis, node)
             else
                dist=0.0
             end if            
               
             call mcD_J_W_F2002(dens,temp,salt,dist)
             call addto(density, node, dens)
          end do
           
       end do
         
       ! reference density is assumed 1 for the pade approximation
       rho_0=1.0

    end if
    
    if (have_option('/material_phase::'//trim(state%name)//'/sediment'))&
         & then 

       call allocate(deltaS, density%mesh, "DeltaS")
       call allocate(remapS, density%mesh, "RemapS")
       call allocate(sedimentdensity, density%mesh, "SedimentDensity")
       ! fluidconcentration is 1-sedimentconcentration.
       call allocate(fluidconcentration, density%mesh, &
            "FluidConcentration")
       call zero(sedimentdensity)
       call set(fluidconcentration, 1.0)
       nSediments=get_nSediments()

       do i=1,nSediments

          class_name = get_sediment_name(i)

          S=>extract_scalar_field(state,trim(class_name))
          option_path = S%option_path

          call get_option(trim(option_path)//'/density', gamma)
          ! Note that 1000.0 is a hack. We actually need to properly
          ! account for the reference density of seawater.
          gamma=rho_0*(gamma/1000.0) - rho_0

          oldS => extract_scalar_field(state, &
               "Old"//trim(class_name))
          
          ! deltaS=theta*S+(1-theta)*oldS-S0
          call remap_field(S, remapS)
          call set(deltaS, remapS)
          call scale(deltaS, theta)
          
          call remap_field(oldS, remapS)
          call addto(deltaS, remapS, 1.0-theta)
          ! density=density+gamma*deltaS
          call addto(sedimentdensity, deltaS, scale=gamma)
          call addto(fluidconcentration, deltaS, scale=-1.0)

       end do
       
       call scale(density,fluidconcentration)
       call addto(density,sedimentdensity)

       call deallocate(deltaS)
       call deallocate(remapS)
       call deallocate(fluidconcentration)
       call deallocate(sedimentdensity)
       
    end if

    if(present(reference_density)) then
      reference_density = rho_0
    end if
    
  end subroutine calculate_perturbation_density

  subroutine mcD_J_W_F2002(density,T,Salinity,distance_to_top)
    !!<  function to evaluate density from the 2002 McDougall, Jackett,
    !!< Wright and Feistel equation of state using Pade approximation.  
    real, intent(out) :: density
    real, intent(in) :: T,Salinity,distance_to_top
   
    real :: p,p1,p2,S

    ! Salinity can be negitive because it's numerically diffused,
    ! some regions may be initialised with zero salinity, and
    ! therefore undershoot may occur.
    S = max(Salinity, 0.0)

    ! calculate pressure in decibars from hydrostatic pressure
    ! using reference density 1000 kg m^-2

    p = 9.81*1000.0*distance_to_top*1.0e-4

    !     evaluate top and bottom of Pade approximant
      
    p1 = 9.99843699e2 &
         + 7.35212840*T - 5.45928211e-2*(T**2) + 3.98476704e-4*(T**3) &
         + 2.96938239*S - 7.23268813e-3*S*T + 2.12382341e-3*(S**2) &
         + 1.04004591e-2*p + 1.03970529e-7*p*(T**2) &
         + 5.18761880e-6*p*S - 3.24041825e-8*(p**2) &
         - 1.23869360e-11*(p**2)*(t**2)

    p2 = 1.0 &
         + 7.28606739e-3*T - 4.60835542e-5*(T**2) + 3.68390573e-7*(T**3) &
         + 1.80809186e-10*(T**4) &
         + 2.14691708e-3*S - 9.27062484e-6*S*T - 1.78343643e-10*S*(T**3) &
         + 4.76534122e-6*(S**1.5) + 1.63410736e-9*(S**1.5)*(T**2) &
         + 5.30848875e-6*p -3.03175128e-16*(p**2)*(t**3) &
         - 1.27934137e-17*(p**3)*T
    
    ! calculate the resulting density
    
    density = p1/p2
    
    ! the perturbation density
    
    density = (density-1000.0)/1000.0
    
  end subroutine mcD_J_W_F2002
  
  subroutine compressible_eos_1mat(state, density, pressure,full_pressure,drhodp,temperature)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional :: density, drhodp
    type(scalar_field), intent(inout), optional :: full_pressure, pressure, temperature

    type(scalar_field), pointer :: hp, complete_pressure, lp
    integer :: stat
    
    character(len=OPTION_PATH_LEN) :: eos_path
    type(scalar_field) :: drhodp_local

    character(len = *), parameter:: hp_name = "HydrostaticPressure"
    character(len = *), parameter:: hpg_name = "HydrostaticPressureGradient"

    ewrite(1,*) 'Entering compressible_eos'
    
    if (present(drhodp)) then
      drhodp_local=drhodp
      if (present(density)) then
         assert(drhodp%mesh==density%mesh)
      end if
      if (present(pressure)) then
         assert(drhodp%mesh==pressure%mesh)
      end if
      if (present(full_pressure)) then
         assert(drhodp%mesh==full_pressure%mesh)
      end if
   else if (present(density)) then
      call allocate(drhodp_local, density%mesh, 'Localdrhop')
   else if (present(pressure)) then
      call allocate(drhodp_local, pressure%mesh, 'Localdrhop')
   else if (present(full_pressure)) then
      call allocate(drhodp_local, full_pressure%mesh, 'Localdrhop')
   else
      FLAbort("No point in being in here if you don't want anything out.")
   end if
   
   eos_path = trim(state%option_path)//'/equation_of_state'
   
   if(have_option(trim(eos_path)//'/compressible')) then
      
      ! each of the following compressible_eos_XXX() routines should always calculate drhodp
      ! (zero if density does not depend on pressure) and calculate density and
      ! pressure if present
      
      if(have_option(trim(eos_path)//'/compressible/stiffened_gas')) then
         
         ! standard stiffened gas eos
         
         if (present(pressure))then
            call compressible_eos_stiffened_gas(state, eos_path,&
                 drhodp_local, &
                 density=density, pressure=pressure)
            else 
               call compressible_eos_stiffened_gas(state, eos_path,&
                 drhodp_local, &
                 density=density, pressure=full_pressure)
            end if
         
      else if(have_option(trim(eos_path)//'/compressible/giraldo2')) then
         
        ! Eq. of state commonly used in atmospheric applications. See
        ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
        ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
        
         if (has_scalar_field(state,hp_name)) then
            hp=>extract_scalar_field(state,hp_name)
            lp=>extract_scalar_field(state,"Pressure")
            allocate(complete_pressure)
            call allocate(complete_pressure,hp%mesh,"FullPressure")
            call remap_field(lp,complete_pressure)
            call addto(complete_pressure,hp)
         else
            complete_pressure=>extract_scalar_field(state,"Pressure")
         end if
         if (present(pressure))then
            call compressible_eos_giraldo_1mat(state, eos_path, drhodp_local, &
                 complete_pressure,density,pressure=pressure,temperature=temperature)
         else
            call compressible_eos_giraldo_1mat(state, eos_path, drhodp_local, &
                 complete_pressure,density=density, pressure=full_pressure,temperature=temperature)
         end if

         if (has_scalar_field(state,hp_name)) then
            call deallocate(complete_pressure)
            deallocate(complete_pressure)
         end if

      elseif(have_option(trim(eos_path)//'/compressible/foam')) then
        
        ! eos used in foam modelling
        
         call compressible_eos_foam(state, eos_path, drhodp_local, &
              density=density, pressure=pressure)

     end if
     
     if (present(pressure) .and. has_scalar_field(state, hp_name)) then
        ! Remove the hydrostatic contribution from full pressure

        hp => extract_scalar_field(state, hp_name)
        call subtract_hydrostatic_pressure_contribution(state,pressure,hp)
     end if
  else
    
     ! I presume we dont' actually want to be here
      FLAbort('Gone into compressible_eos without having equation_of_state/compressible')


    end if

    if(present(density)) then
      ewrite_minmax(density)
    end if

    if(present(pressure)) then
      ewrite_minmax(pressure)
    end if

    if(present(full_pressure)) then
      ewrite_minmax(full_pressure%val)
    end if

    if(present(drhodp)) then      
      ewrite_minmax(drhodp)
    else
      call deallocate(drhodp_local)
    end if

  end subroutine compressible_eos_1mat

  subroutine compressible_eos_mmat(state, density, pressure,full_pressure,drhodp,temperature)

    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), intent(inout), optional :: density, drhodp
    type(scalar_field), intent(inout), optional :: full_pressure, pressure, temperature

    type(scalar_field), pointer :: hp, complete_pressure, lp
    integer :: stat
    
    character(len=OPTION_PATH_LEN) :: eos_path
    type(scalar_field) :: drhodp_local

    if (size(state)==1) then
       call compressible_eos_1mat(state(1), density, pressure,full_pressure,drhodp,temperature)
       return
    end if

    ewrite(1,*) 'Entering compressible_eos'
    
    if (present(drhodp)) then
      drhodp_local=drhodp
      if (present(density)) then
         assert(drhodp%mesh==density%mesh)
      end if
      if (present(pressure)) then
         assert(drhodp%mesh==pressure%mesh)
      end if
      if (present(full_pressure)) then
         assert(drhodp%mesh==full_pressure%mesh)
      end if
   else if (present(density)) then
      call allocate(drhodp_local, density%mesh, 'Localdrhop')
   else if (present(pressure)) then
      call allocate(drhodp_local, pressure%mesh, 'Localdrhop')
   else if (present(full_pressure)) then
      call allocate(drhodp_local, full_pressure%mesh, 'Localdrhop')
   else
      FLAbort("No point in being in here if you don't want anything out.")
   end if
   
   eos_path = trim(state(1)%option_path)//'/equation_of_state'
   
   if(have_option(trim(eos_path)//'/compressible')) then       
      if(have_option(trim(eos_path)//'/compressible/ATHAM')) then
         
        ! Eq. of state commonly used in atmospheric applications. See
        ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
        ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
        
         complete_pressure=>extract_scalar_field(state(1),"Pressure")
         if (present(pressure))then
            call compressible_eos_atham(state, eos_path, drhodp_local, &
                 complete_pressure,density,pressure=pressure,temperature=temperature)
         else
            call compressible_eos_atham(state, eos_path, drhodp_local, &
                 complete_pressure,density=density, pressure=full_pressure,temperature=temperature)
         end if

      else

               FLAbort('Gone into multimaterial compressible_eos without having equation_of_state/compressible/ATHAM')

      end if
     
  else
    
     ! I presume we dont' actually want to be here
      FLAbort('Gone into compressible_eos without having equation_of_state/compressible')


    end if

    if(present(density)) then
      ewrite_minmax(density)
    end if

    if(present(pressure)) then
      ewrite_minmax(pressure)
    end if

    if(present(full_pressure)) then
      ewrite_minmax(full_pressure%val)
    end if

    if(present(drhodp)) then      
      ewrite_minmax(drhodp)
    else
      call deallocate(drhodp_local)
    end if

  end subroutine compressible_eos_mmat
    
  subroutine compressible_eos_stiffened_gas(state, eos_path, drhodp, &
    density, pressure)
    ! Standard stiffened gas equation
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure
    
    !locals
    integer :: stat, gstat, cstat
    type(scalar_field), pointer :: pressure_local, energy_local, density_local, hp
    real :: reference_density, ratio_specific_heats
    real :: bulk_sound_speed_squared, atmospheric_pressure
    type(scalar_field) :: energy_remap, pressure_remap, density_remap
    logical :: incompressible

    character(len = *), parameter:: hp_name = "HydrostaticPressure"
    character(len = *), parameter:: hpg_name = "HydrostaticPressureGradient"
    
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/reference_density', &
                        reference_density, default=0.0)
        
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/ratio_specific_heats', &
                    ratio_specific_heats, stat=gstat)
    if(gstat/=0) then
      ratio_specific_heats=1.0
    end if
    
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/bulk_sound_speed_squared', &
                    bulk_sound_speed_squared, stat=cstat)
    if(cstat/=0) then
      bulk_sound_speed_squared=0.0
    end if
    
    incompressible = ((gstat/=0).and.(cstat/=0))
    if(incompressible) then
      ewrite(0,*) "Selected compressible eos but not specified a bulk_sound_speed_squared or a ratio_specific_heats."
    end if
    
    call zero(drhodp)
    
    if(.not.incompressible) then
      energy_local=>extract_scalar_field(state,'InternalEnergy',stat=stat)
      ! drhodp = 1.0/( bulk_sound_speed_squared + (ratio_specific_heats - 1.0)*energy )
      if((stat==0).and.(gstat==0)) then   ! we have an internal energy field and we want to use it
        call allocate(energy_remap, drhodp%mesh, 'RemappedInternalEnergy')
        call remap_field(energy_local, energy_remap)
        
        call addto(drhodp, energy_remap, (ratio_specific_heats-1.0))
        
        call deallocate(energy_remap)
      end if
      call addto(drhodp, bulk_sound_speed_squared)
      call invert(drhodp)
    end if

    if(present(density)) then
      ! calculate the density
      ! density may equal density in state depending on how this
      ! subroutine is called
      if(incompressible) then
        ! density = reference_density
        call set(density, reference_density)
      else
        pressure_local=>extract_scalar_field(state,'Pressure',stat=stat)
        if (has_scalar_field(state,hp_name)) then
           hp => extract_scalar_field(state,hp_name)
           call addto(pressure_local,hp)
        end if
        if (stat==0) then
          assert(density%mesh==drhodp%mesh)
        
          ! density = drhodp*(pressure_local + atmospheric_pressure
          !                  + bulk_sound_speed_squared*reference_density)
          call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                          atmospheric_pressure, default=0.0)
          
          call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
          call remap_field(pressure_local, pressure_remap)
          
          call set(density, reference_density*bulk_sound_speed_squared + atmospheric_pressure)
          call addto(density, pressure_remap)
          call scale(density, drhodp)
          
          call deallocate(pressure_remap)
        else
          FLExit('No Pressure in material_phase::'//trim(state%name))
        end if
      end if
    end if

    if(present(pressure)) then
      if(incompressible) then
        ! pressure is unrelated to density in this case
        call zero(pressure)
      else
        ! calculate the pressure using the eos and the calculated (probably prognostic)
        ! density
        density_local=>extract_scalar_field(state,'Density',stat=stat)
        if (stat==0) then
          assert(pressure%mesh==drhodp%mesh)
          
          ! pressure = density_local/drhodp &
          !          - bulk_sound_speed_squared*reference_density
          
          call allocate(density_remap, drhodp%mesh, "RemappedDensity")
          call remap_field(density_local, density_remap)
          
          call set(pressure, drhodp)
          call invert(pressure)
          call scale(pressure, density_remap)
          call addto(pressure, -bulk_sound_speed_squared*reference_density)
          
          call deallocate(density_remap)
        else
          FLExit('No Density in material_phase::'//trim(state%name))
        end if
      end if
    end if

  end subroutine compressible_eos_stiffened_gas
  
  subroutine compressible_eos_giraldo(state, eos_path, drhodp, &
    complete_pressure,density, pressure)
    ! Eq. of state commonly used in atmospheric applications. See
    ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
    ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure
    type(scalar_field), intent(in) :: complete_pressure
    
    ! locals
    integer :: stat, gstat, cstat, pstat, tstat
    type(scalar_field), pointer :: pressure_local, energy_local, density_local, &
                                   & temperature_local
    type(scalar_field) :: energy_remap, pressure_remap, density_remap, &
                          & temperature_remap
    real :: reference_density, p_0, c_p, c_v
    real :: drhodp_node, power
    real :: R
!     type(scalar_field) :: pressure_remap, density_remap, temperature_remap
    logical :: incompressible
    integer :: node
    
    call get_option(trim(eos_path)//'/compressible/giraldo2/reference_pressure', &
                    p_0, default=1.0e5)
        
    call get_option(trim(eos_path)//'/compressible/giraldo2/C_P', &
                    c_p, stat=gstat)
    if(gstat/=0) then
       c_p=1000.0
    end if
        
    call get_option(trim(eos_path)//'/compressible/giraldo2/C_V', &
                    c_v, stat=cstat)
    if(cstat/=0) then
       c_v=714.285714
    end if

!     print *, 'In giraldo2'
!     p_0=1.0e5
!     c_p=1000.0
!     c_v=714.285714

!     call get_option(trim(eos_path)//'/compressible/giraldo2/reference_height', &
!                     ref_h, stat=hstat)
!     if(hstat/=0) then
!        ref_h=0.0
!     end if

    R=c_p-c_v
        
    incompressible = ((gstat/=0).or.(cstat/=0))
    if(incompressible) then
       ewrite(0,*) "Selected compressible eos but not specified either C_P or C_V."
    end if

!     if (hstat/=0) then
!         ewrite(0,*) "Warning: No reference height set for the equation of state. Defaulting to 0."
!     end if

        ! Were going to need the position and velocity fields for various calculations.
        ! Currently some dirty, dirty remaps will be going on, but as this is only being
        ! tested with linear shape functions for all spaces it's ok for now.

!     x_local => extract_vector_field(state, "Coordinate")
!         u_local => extract_vector_field(state, "NonlinearVelocity")
!         u_local => extract_vector_field(state, "Velocity")

!         call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, stat)

!         dim=x_local%dim
!         assert(x_local%dim==u_local%dim)
! 
!         allocate(x(dim), u(dim))
        
        call zero(drhodp)
        
        if(.not.incompressible) then
          energy_local=>extract_scalar_field(state,'InternalEnergy',stat=stat)

          if((gstat==0).and.(cstat==0)) then
            
            call allocate(energy_remap, drhodp%mesh, 'RemappedInternalEnergy')
!             call allocate(x_remap, x_local%dim, drhodp_local%mesh, name="RemappedCoordinate")
!             call allocate(u_remap, u_local%dim, drhodp_local%mesh, name="RemappedVelocity")

            call remap_field(energy_local, energy_remap)
!             x_remap = get_nodal_coordinate_field( state, drhodp_local%mesh )
!             call remap_field(u_local, u_remap, stat=stat)
            
            do node=1,node_count(drhodp)
!               x=node_val(x_remap,node)
!               u=node_val(u_remap,node)
!               drhodp_node=(R/c_v)*(node_val(energy_remap,node)-0.5*dot_product(u,u)-gravity_magnitude*(x(dim)-ref_h))
              drhodp_node=(R/c_v)*node_val(energy_remap,node)
              call set(drhodp, node, drhodp_node)
            end do
            
            call deallocate(energy_remap)
!             call deallocate(x_remap)
!             call deallocate(u_remap)
          end if
          call invert(drhodp)

        end if

        if(present(density)) then
          ! calculate the density
          ! density may equal density in state depending on how this
          ! subroutine is called
          if(incompressible) then
            ! density = reference_density
            call set(density, reference_density)
          else
              call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
              call remap_field(complete_pressure, pressure_remap)

!               call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
!                               atmospheric_pressure, default=0.0)
              
!               call set(density, atmospheric_pressure)
!               call addto(density, pressure_remap)
              call set(density, pressure_remap)
              call scale(density, drhodp)
              
              call deallocate(pressure_remap)
           end if
        end if

       if(present(pressure)) then
          if(incompressible) then
            ! pressure is unrelated to density in this case
            call zero(pressure)
          else
            ! calculate the pressure using the eos and the calculated (probably prognostic)
            ! density
             
            density_local=>extract_scalar_field(state,'Density',stat=stat)

            if (stat==0) then
              assert(pressure%mesh==drhodp%mesh)
              
              call allocate(density_remap, drhodp%mesh, "RemappedDensity")
              call remap_field(density_local, density_remap)
              
              call set(pressure, drhodp)
              call invert(pressure)
              call scale(pressure, density_remap)

              call deallocate(density_remap)
            else
              FLExit('No Density in material_phase::'//trim(state%name))
            end if
          end if
        end if

!         deallocate(x, u)

  end subroutine compressible_eos_giraldo

  subroutine make_bulk_quantities(bulk,q,q_val,denom)
    type(scalar_field), intent(inout) :: bulk
    type(scalar_field), intent(in) :: q(:)
    type(scalar_field), intent(inout), optional :: denom
    real, intent(in) :: q_val(size(q))
    type(scalar_field) :: fac
    
    integer :: i
    
    ! Warning, this routine is a temporary hack to get bulk quantities
    
    
    
    call zero(bulk)
    do i=1,size(q)
       call addto(bulk,q(i),scale=q_val(i))
    end do
    
    if (present(denom)) then
       call allocate(fac,denom%mesh,"scale_factor")
       call set(fac,denom)
       call invert(fac)
       call scale(bulk,fac)
       call deallocate(fac)
    end if
    
    
  end subroutine make_bulk_quantities
  
  subroutine make_bulk_quantities_mmat(bulk,state,vtype)
    type(scalar_field), intent(inout) :: bulk
    type(state_type), intent(in), dimension(:) :: state(:)
    character(len=*), intent(in) :: vtype
    type(scalar_field), pointer :: fraction
    type(scalar_field) :: dry_fraction
    
    integer :: i, stat
    real :: q_val
    character(len=OPTION_PATH_LEN) :: eos_path
    
    ! Warning, this routine is a temporary hack to get bulk quantities
    
    call zero(bulk)
    call allocate(dry_fraction,bulk%mesh,"DryGasMassFraction")
    call set(dry_fraction,1.0)
    do i=2,size(state)
       if (has_scalar_field(state(i),"MassFraction")) then
          eos_path = trim(state(i)%option_path)//'/equation_of_state'
          call get_option(trim(eos_path)//'/compressible/ATHAM/'//vtype, &
               q_val, stat=stat)
          fraction=>extract_scalar_field(state(i),"MassFraction")
          if (stat /= 0) then 
             ewrite(0,*) "Selected compressible eos but not specified "//vtype//" in "//state(i)%name//"."
          else
             call addto(bulk,fraction,scale=q_val)
             call addto(dry_fraction,fraction,scale=-1.0)
          end if
       end if
    end do
    eos_path = trim(state(1)%option_path)//'/equation_of_state'
    call get_option(trim(eos_path)//'/compressible/ATHAM/'//vtype, &
         q_val, stat=stat)
    if (stat /= 0) then 
       ewrite(0,*) "Selected compressible eos but not specified "//vtype//" in "//state(1)%name//"."
    else
       call addto(bulk,dry_fraction,scale=q_val)
    end if
    call deallocate(dry_fraction)

    print*, vtype, minval(bulk%val), maxval(bulk%val)
    
  end subroutine make_bulk_quantities_mmat

  subroutine make_incompressible_fix_mmat(incompfix,state,density)
    type(scalar_field), intent(inout) :: incompfix
    type(state_type), intent(in), dimension(:) :: state
    type(scalar_field), intent(in) :: density

    type(scalar_field), pointer :: fraction
    character(len=OPTION_PATH_LEN) :: eos_path
    integer :: i,stat
    real :: rho

    call zero(incompfix)

    do i=1,size(state)
       if (has_scalar_field(state(i),"MassFraction")) then
          eos_path = trim(state(i)%option_path)//'/equation_of_state'
          call get_option(trim(eos_path)//'/compressible/ATHAM/density', &
               rho, stat=stat)
          if (stat == 0) then
             fraction=>extract_scalar_field(state(i),"MassFraction")
             call addto(incompfix,fraction,scale=-1.0/rho)
          end if
       end if
    end do
    call scale(incompfix,density)
    call addto(incompfix,1.0)

  end subroutine make_incompressible_fix_mmat

  subroutine make_atham_quantities(state,qc_p,qc_v,qc_solid,scaledq)

    type(state_type), intent(in), dimension(:) :: state
    type(scalar_field), intent(inout) :: qc_p,qc_v,qc_solid,scaledq

    type(scalar_field), pointer :: fraction
    character(len=OPTION_PATH_LEN) :: eos_path
    integer :: i,stat
    real :: c_v,c_p,rho
    type(scalar_field) :: dry_fraction

    call zero(qc_p)
    call zero(qc_v)
    call zero(qc_solid)
    call zero(scaledq)

    call allocate(dry_fraction,qc_solid%mesh,"DryGasMassFraction")
    call set(dry_fraction,1.0)

    do i=2,size(state)
       if (has_scalar_field(state(i),"MassFraction")) then
          fraction=>extract_scalar_field(state(i),"MassFraction")
          call addto(dry_fraction,fraction,scale=-1.0)
          eos_path = trim(state(i)%option_path)//'/equation_of_state'
          call get_option(trim(eos_path)//'/compressible/ATHAM/C_P', &
               c_p, stat=stat)
          if (stat /= 0) then
             ewrite(0,*) "Selected compressible eos but not specified C_P in "//state(i)%name//"."
             cycle
          end if
          call get_option(trim(eos_path)//'/compressible/ATHAM/C_V', &
               c_v, stat=stat)
          if (stat /= 0) then 
             ewrite(0,*) "Selected compressible eos but not specified C_V in "//state(i)%name//"."
             cycle
          end if
          if (have_option(trim(eos_path)//'/compressible/ATHAM/density')) then
             call get_option(trim(eos_path)//'/compressible/ATHAM/density',rho,stat)
             call addto(qc_solid,fraction,scale=c_v)
             call addto(scaledq,fraction,scale=1.0/rho)
          else
             call addto(qc_v,fraction,scale=c_v)
             call addto(qc_p,fraction,scale=c_p)
          end if
       end if
    end do
    eos_path = trim(state(1)%option_path)//'/equation_of_state'
    call get_option(trim(eos_path)//'/compressible/ATHAM/C_P', &
         c_p, stat=stat)
    if (stat /= 0) then 
       ewrite(0,*) "Selected compressible eos but not specified C_P in "//state(i)%name//"."
       c_p=0.0
    end if
    call get_option(trim(eos_path)//'/compressible/ATHAM/C_V', &
         c_v, stat=stat)
    if (stat /= 0) then 
       ewrite(0,*) "Selected compressible eos but not specified C_V in "//state(i)%name//"."
       c_v=0.0
    end if
    call addto(qc_v,dry_fraction,scale=c_v)
    call addto(qc_p,dry_fraction,scale=c_p)

    call deallocate(dry_fraction)
    
  end subroutine make_atham_quantities


  subroutine compressible_eos_atham(state, eos_path, drhodp, &
    complete_pressure,density, pressure,temperature)
    ! Eq. of state commonly used in atmospheric applications. See
    ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
    ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
    type(state_type), intent(inout), dimension(:) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure,&
         temperature
    type(scalar_field), intent(in) :: complete_pressure
    
    ! locals
    integer :: stat, gstat, cstat, pstat, tstat
    type(scalar_field), pointer :: pressure_local, energy_local,&
         & density_local, &
         & temperature_local
    type(scalar_field) :: energy_remap, pressure_remap, density_remap, &
                          & temperature_remap,qc_v,qc_p,qc_solid,scaledq,&
                          q_g, incompfix, temp_local, rhs, dpress
    type(vector_field), pointer :: X
    real :: reference_density, p0, c_p, c_v
    real :: drhodp_node, power, temperature_node, density_node,&
         pressure_node, energy_node,&
         c_v_node, R_node, q_node, d_node, i_node
    real :: R, c_p_water, c_v_water, R_water, rho_w
    logical :: incompressible
    integer :: node, thermal_variable, ele, nlin
    type(csr_matrix) :: M
    type(csr_sparsity), pointer :: M_sparsity
    
    call get_option(trim(eos_path)//'/compressible/ATHAM/reference_pressure', &
                    p0, default=1.0e5)

    X=>extract_vector_field(state(1), "Coordinate")

    call allocate(qc_p,drhodp%mesh,"qc_p")
    call allocate(qc_v,drhodp%mesh,"qc_v")
    call allocate(qc_solid,drhodp%mesh,"qc_solid")
    call allocate(scaledq,drhodp%mesh,"scaledp")

    if (present(temperature)) then
       call allocate(temp_local,temperature%mesh,"Temperature")
    end if

    call make_atham_quantities(state,qc_p,qc_v,qc_solid,scaledq)


    density_local=>extract_scalar_field(state(1),"Density",stat=cstat)
    pressure_local=>extract_scalar_field(state(1),"Pressure",stat=cstat)
    call zero(drhodp)

    if (drhodp%mesh==density_local%mesh) then 
       if(drhodp%mesh==pressure_local%mesh ) then
          call set(drhodp,pressure_local)
       else
          call calculate_galerkin_projection(state(1),pressure_local,drhodp)
       end if
       call invert(drhodp)
       call scale(drhodp,density_local)
    elseif (drhodp%mesh==pressure_local%mesh) then
       call calculate_galerkin_projection(state(1),density_local,drhodp)
       call invert(drhodp)
       call scale(drhodp,pressure_local)
       call invert(drhodp)
    else
       FLAbort("drhodp not on pressure or density mesh!")
    end if
    
    if (has_scalar_field(state(1),'InternalEnergy')) then
       thermal_variable=1
       energy_local=>extract_scalar_field(state(1),'InternalEnergy',stat=stat)
    else if (has_scalar_field(state(1),'PotentialTemperature')) then
       thermal_variable=2
       energy_local=>extract_scalar_field(state(1),'PotentialTemperature',&
            stat=stat)
       call allocate(dpress,drhodp%mesh,"LocalDeltaP")
    else
       FLAbort("No thermodynamic variable!")
    end if

    M_sparsity=>get_csr_sparsity_firstorder(state(1),drhodp%mesh,drhodp%mesh)
    call allocate(M,M_sparsity,name="EOSMatrix")
    call allocate(rhs,drhodp%mesh,name="RHS")

    call zero(M)
    call zero(rhs)


    do ele=1,ele_count(drhodp)
       call assemble_drhodp_matrix(ele,M,rhs,X,drhodp,&
            density_local,pressure_local,&
            energy_local,qc_v,qc_p,qc_solid,scaledq,p0,thermal_variable)
    end do


    ewrite_minmax(rhs)
    ewrite_minmax(drhodp)
    ewrite_minmax(density_local)
    ewrite_minmax(pressure_local)
    drhodp%option_path=pressure_local%option_path
    call petsc_solve(drhodp,M,rhs)
          

    if(present(density)) then
       ! calculate the density
       ! density may equal density in state depending on how this
       ! subroutine is called

       call zero(M)
       call zero(rhs)
!       call remap_field(density_local,density)
       call zero(density)

       do ele=1,ele_count(density)
          call assemble_density_matrix(ele,M,rhs,X,density,&
               pressure_local,&
               energy_local,qc_v,qc_p,qc_solid,scaledq,p0,thermal_variable)
       end do

       call petsc_solve(density,M,rhs,&
            option_path=pressure_local%option_path)
    end if
    
    if(present(pressure)) then
       ! calculate the pressure using the eos and the calculated (probably prognostic)
       ! density
       
       if (stat==0) then
          assert(pressure%mesh==drhodp%mesh)
          
!          call remap_field(complete_pressure,pressure)
!          call set(pressure,p0)
          call calculate_galerkin_projection(state(1),pressure_local,pressure)

          select case(thermal_variable)
          case(1)
             call zero(M)
             call zero(rhs)


             do ele=1,ele_count(pressure)
                call assemble_pressure_matrix(ele,M,rhs,X,density_local,&
                     pressure,&
                     energy_local,qc_v,qc_p,qc_solid,&
                     scaledq,p0,thermal_variable)
             end do


             call petsc_solve(pressure,M,rhs,&
                  option_path=pressure_local%option_path)

          case(2)

             if (minval(pressure%val)<1e-7)&
                  call set(pressure,p0)

             do nlin=1,10

                call zero(M)
                call zero(rhs)
                call zero(dpress)

                do ele=1,ele_count(pressure)
                   call assemble_pressure_matrix(ele,M,rhs,X,density_local,&
                        pressure,energy_local,qc_v,qc_p,qc_solid,&
                        scaledq,p0,thermal_variable)
                end do

                call petsc_solve(dpress,M,rhs,&
                     option_path=pressure_local%option_path)

                call addto(pressure,dpress)

             end do
          end select
       else
          FLExit('No Density in material_phase::'//trim(state(1)%name))
       end if
    end if
    

        if (present(temperature)) then
       call zero(M)
       call zero(rhs)
       call remap_field(temperature,temp_local)

       do ele=1,ele_count(temp_local)
          call assemble_temperature_matrix(ele,M,rhs,X,temp_local,&
               density_local,complete_pressure,&
               energy_local,qc_v,qc_p,qc_solid,scaledq,p0,thermal_variable)
       end do

    
       temp_local%option_path=pressure_local%option_path
       call petsc_solve(temp_local,M,rhs)

    end if


    call deallocate(qc_p)
    call deallocate(qc_v)
    call deallocate(qc_solid)
    call deallocate(scaledq)
    call deallocate(M)
    call deallocate(rhs)
    if (present(temperature)) then
       call set(temperature,temp_local)
       call deallocate(temp_local)
    end if
    if (thermal_variable==2) call deallocate(dpress)

  end subroutine compressible_eos_atham
  
  subroutine assemble_drhodp_matrix(ele,M,rhs,x,drhodp,&
       density,pressure,energy,&
       qc_v,qc_p,qc_solid,scaledq,p0,thermal_variable)
    
    integer, intent(in) :: ele, thermal_variable
    real, intent(in) :: p0
    type(csr_matrix), intent(inout) :: M
    type(scalar_field) :: rhs
    type(vector_field), intent(in) :: X
    type(scalar_field), intent(in) :: drhodp,density,pressure,energy,&
         qc_v,qc_p,qc_solid,scaledq
    type(element_type), pointer :: myshape
    integer, dimension(:), pointer :: nodes
    real, dimension(ele_ngi(drhodp,ele)) :: detwei

    real, dimension(ele_ngi(drhodp, ele)) :: qcv_gi,&
         qcp_gi, qcs_gi, R_gi, p_gi
    

    
    nodes=> ele_nodes(drhodp,ele)
    myshape=> ele_shape(drhodp,ele)

    qcv_gi=ele_val_at_quad(qc_v,ele)
    qcp_gi=ele_val_at_quad(qc_p,ele)
    qcs_gi=ele_val_at_quad(qc_solid,ele)
    R_gi=qcp_gi-qcv_gi
    
    call transform_to_physical(X, ele,detwei)
    
    select case(thermal_variable)

    case(1)
       
       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
            R_gi*ele_val_at_quad(energy,ele)*detwei))
       call addto(rhs,nodes,shape_rhs(myshape,&
            (1.0-ele_val_at_quad(density,ele)&
            *ele_val_at_quad(scaledq,ele))**2&
            *(qcv_gi+qcs_gi)*detwei))
    case(2)
       
       p_gi=ele_val_at_quad(pressure,ele)

       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
            p_gi**(R_gi/qcp_gi)*(R_gi*(qcp_gi+qcs_gi)*ele_val_at_quad(energy,ele)&
            +ele_val_at_quad(scaledq,ele)&
            *(p0**(R_gi/qcp_gi)*(qcp_gi)+p_gi**(R_gi/qcp_gi)*qcs_gi))&
            *detwei))
       call addto(rhs,nodes,shape_rhs(myshape,&
            (1.0-ele_val_at_quad(density,ele)&
            *ele_val_at_quad(scaledq,ele))&
            *(p0**(R_gi/qcp_gi)*(qcv_gi)+p_gi**(R_gi/qcp_gi)*qcs_gi)&
            *detwei))

    end select
       
    
  end subroutine assemble_drhodp_matrix
  
  subroutine assemble_temperature_matrix(ele,M,rhs,x,temp,&
       density,pressure,energy,&
       qc_v,qc_p,qc_solid,scaledq,p0,thermal_variable)
    
    integer, intent(in) :: ele, thermal_variable
    real, intent(in) :: p0
    type(csr_matrix), intent(inout) :: M
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: X
    type(scalar_field), intent(in) :: temp,density,pressure,energy,&
         qc_v,qc_p,qc_solid,scaledq
    type(element_type), pointer :: myshape
    integer, dimension(:), pointer :: nodes
    real, dimension(ele_ngi(temp,ele)) :: detwei

    real, dimension(ele_ngi(temp, ele)) :: qcv_gi,&
         qcp_gi, qcs_gi, R_gi, p_gi
    
    
    nodes=> ele_nodes(temp,ele)
    myshape=> ele_shape(temp,ele)

    qcv_gi=ele_val_at_quad(qc_v,ele)
    qcp_gi=ele_val_at_quad(qc_p,ele)
    qcs_gi=ele_val_at_quad(qc_solid,ele)
    R_gi=qcp_gi-qcv_gi
    
    call transform_to_physical(X, ele, detwei)
    
    select case(thermal_variable)

    case(1)

       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
            (ele_val_at_quad(qc_v,ele)+ele_val_at_quad(qc_solid,ele))*detwei))
       call addto(rhs,nodes,shape_rhs(myshape,&
            ele_val_at_quad(energy,ele)*detwei))

    case(2)

       p_gi=ele_val_at_quad(pressure,ele)
       
       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
            (p0**(R_gi/qcp_gi)*qcp_gi+p_gi**(R_gi/qcp_gi)*qcs_gi)*detwei))
       call addto(rhs,nodes,shape_rhs(myshape,&
            p_gi**(R_gi/qcp_gi)*(qcp_gi+qcs_gi)&
            *ele_val_at_quad(energy,ele)*detwei))

    end select
       
        
  end subroutine assemble_temperature_matrix

subroutine assemble_density_matrix(ele,M,rhs,x,&
       density,pressure,energy,&
       qc_v,qc_p,qc_solid,scaledq,p0,thermal_variable)
    
    integer, intent(in) :: ele, thermal_variable
    type(csr_matrix), intent(inout) :: M
    real, intent(in) :: p0
    type(scalar_field) :: rhs
    type(vector_field), intent(in) :: X
    type(scalar_field), intent(in) :: density,pressure,energy,&
         qc_v,qc_p,qc_solid,scaledq
    type(element_type), pointer :: myshape
    integer, dimension(:), pointer :: nodes
    real, dimension(ele_ngi(density,ele)) :: detwei

    real, dimension(ele_ngi(density, ele)) :: qcv_gi,&
         qcp_gi, qcs_gi, R_gi, p_gi
    
    
    nodes=> ele_nodes(density,ele)
    myshape=> ele_shape(density,ele)

    qcv_gi=ele_val_at_quad(qc_v,ele)
    qcp_gi=ele_val_at_quad(qc_p,ele)
    qcs_gi=ele_val_at_quad(qc_solid,ele)
    R_gi=qcp_gi-qcv_gi
    p_gi=ele_val_at_quad(pressure,ele)
    
    call transform_to_physical(X, ele,detwei)
    
    select case(thermal_variable)

    case(1)

       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
            (R_gi&
            *ele_val_at_quad(energy,ele)&
            +(qcv_gi+qcs_gi)*&
            p_gi*ele_val_at_quad(scaledq,ele))*detwei))
       call addto(rhs,nodes,shape_rhs(myshape,&
            (qcv_gi+qcs_gi)*p_gi*detwei))

    case(2)

       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
            (R_gi*(qcp_gi+qcs_gi)*ele_val_at_quad(energy,ele)+&
            p_gi*((p0/p_gi)**(R_gi/qcp_gi)*qcp_gi&
            +qcs_gi)*ele_val_at_quad(scaledq,ele))*detwei))
       call addto(rhs,nodes,shape_rhs(myshape,&
            p_gi*((p0/p_gi)**(R_gi/qcp_gi)*qcp_gi+qcs_gi)*detwei))
       
    
    end select

  end subroutine assemble_density_matrix

subroutine assemble_pressure_matrix(ele,M,rhs,x,&
       density,pressure,energy,&
       qc_v,qc_p,qc_solid,scaledq,p0,thermal_variable)
    
    integer, intent(in) :: ele, thermal_variable
    real, intent(in) :: p0
    type(csr_matrix), intent(inout) :: M
    type(scalar_field) :: rhs
    type(vector_field), intent(in) :: X
    type(scalar_field), intent(in) :: density,pressure,energy,&
         qc_v,qc_p,qc_solid,scaledq
    type(element_type), pointer :: myshape
    integer, dimension(:), pointer :: nodes
    real, dimension(ele_ngi(pressure,ele)) :: detwei
    
    real, dimension(ele_ngi(pressure, ele)) :: qcv_gi,&
         qcp_gi, qcs_gi, R_gi, p_gi, rho_gi, dqs_gi
    
    nodes=> ele_nodes(pressure,ele)
    myshape=> ele_shape(pressure,ele)

    qcv_gi=ele_val_at_quad(qc_v,ele)
    qcp_gi=ele_val_at_quad(qc_p,ele)
    qcs_gi=ele_val_at_quad(qc_solid,ele)
    R_gi=qcp_gi-qcv_gi
    p_gi=ele_val_at_quad(pressure,ele)
    rho_gi=ele_val_at_quad(density,ele)
    dqs_gi=ele_val_at_quad(scaledq,ele)
    
    call transform_to_physical(X, ele,detwei)
    
    select case(thermal_variable)

    case(1)

       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
            (qcv_gi+qcs_gi)*(1.0-rho_gi*dqs_gi)*detwei))
       call addto(rhs,nodes,shape_rhs(myshape,&
            R_gi*rho_gi*ele_val_at_quad(energy,ele)*detwei))

    case(2)

!       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
!            ((1.0-rho_gi*dqs_gi)*p_gi**(qcv_gi/qcp_gi)*&
!            (p0**(R_gi/qcp_gi)*(qcp_gi+qcv_gi)&
!            +2.0*p_gi**(R_gi/qcp_gi)*qcs_gi)&
!            -rho_gi*R_gi*(qcp_gi+qcs_gi)*ele_val_at_quad(energy,ele))*detwei))

       call addto(M,nodes,nodes,shape_shape(myshape,myshape,&
            (1.0-rho_gi*dqs_gi)*((p0/p_gi)**(R_gi/qcp_gi)*qcv_gi+qcs_gi)&
            *detwei))
!       call addto(rhs,nodes,shape_rhs(myshape,&
!            p_gi*(rho_gi*R_gi*(qcp_gi+qcs_gi)*ele_val_at_quad(energy,ele)&
!            -(1.0-rho_gi*dqs_gi)*p_gi**(qcv_gi/qcp_gi)&
!            *(p0**(R_gi/qcp_gi)*qcp_gi+p_gi**(R_gi/qcp_gi)*qcs_gi))*detwei))
       call addto(rhs,nodes,shape_rhs(myshape,&
            (rho_gi*R_gi*(qcp_gi+qcs_gi)*ele_val_at_quad(energy,ele)&
            -(1.0-rho_gi*dqs_gi)*p_gi&
            *((p0/p_gi)**(R_gi/qcp_gi)*qcp_gi+qcs_gi))*detwei))


    end select
    
  end subroutine assemble_pressure_matrix
      
  subroutine compressible_eos_giraldo_1mat(state, eos_path, drhodp, &
    complete_pressure,density, pressure,temperature)
    ! Eq. of state commonly used in atmospheric applications. See
    ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
    ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure,&
         temperature
    type(scalar_field), intent(in) :: complete_pressure
    
    ! locals
    integer :: stat, gstat, cstat, pstat, tstat
    type(scalar_field), pointer :: pressure_local, energy_local,&
         & density_local, &
         & temperature_local, q_v,q_c,q_r
    type(scalar_field) :: energy_remap, pressure_remap, density_remap, &
                          & temperature_remap, R_bulk,c_p_bulk,c_v_bulk,&
                          q_g, incompfix, temp_local
    type(scalar_field), target :: dummyscalar
    real :: reference_density, p_0, c_p, c_v
    real :: drhodp_node, power, temperature_node, density_node,&
         pressure_node, energy_node,&
         c_v_node, R_node, q_node, d_node, i_node
    real :: R, c_p_water, c_v_water, R_water, rho_w
    logical :: incompressible
    integer :: node

    character(len = *), parameter:: hp_name = "HydrostaticPressure"
    character(len = *), parameter:: hpg_name = "HydrostaticPressureGradient"
    
    call get_option(trim(eos_path)//'/compressible/giraldo2/reference_pressure', &
                    p_0, default=1.0e5)
        
    call get_option(trim(eos_path)//'/compressible/giraldo2/C_P', &
                    c_p, stat=gstat)
    if(gstat/=0) then
       c_p=1000.0
    end if
        
    call get_option(trim(eos_path)//'/compressible/giraldo2/C_V', &
                    c_v, stat=cstat)
    if(cstat/=0) then
       c_v=714.285714
    end if

    R=c_p-c_v

    c_p_water=1840
    c_v_water=1840-461.5
    R_water=c_p_water-c_v_water
    rho_w=1000.0

    call allocate(dummyscalar,drhodp%mesh,"dummy")
    call zero(dummyscalar)

    call allocate(R_bulk,drhodp%mesh,"BulkR")
    call allocate(c_v_bulk,drhodp%mesh,"Bulkc_v")

    if (has_scalar_field(state,"WaterVapourFraction")) then
       q_v=>extract_scalar_field(state,"WaterVapourFraction")
    else
       q_v=>dummyscalar
    end if

    if (has_scalar_field(state,"CloudWaterFraction")) then
       q_c=>extract_scalar_field(state,"CloudWaterFraction")
    else
       q_c=>dummyscalar
    end if

    if (has_scalar_field(state,"RainWaterFraction")) then
       q_r=>extract_scalar_field(state,"RainWaterFraction")
    else
       q_r=>dummyscalar
    end if

    if (present(temperature)) then
       call allocate(temp_local,q_r%mesh,"Temperature")
    end if

    call allocate(q_g,drhodp%mesh,"GasMassFraction")
    call allocate(incompfix,drhodp%mesh,"IncompressibleScaleFactor2")

    call set (q_g,1.0)
    call addto(q_g,q_v,scale=-1.0)
    call addto(q_g,q_c,scale=-1.0)
    call addto(q_g,q_r,scale=-1.0)


    call make_bulk_quantities(R_bulk,(/q_g,q_v/),&
         (/R,R_water/))
    call make_bulk_quantities(c_v_bulk,(/q_g,q_v,q_c,q_r/),&
         (/c_v,c_v_water,4100.0,4100.0/))

    density_local=>extract_scalar_field(state,"Density",stat=cstat)
    pressure_local=>extract_scalar_field(state,"Pressure",stat=cstat)

    call zero(incompfix)
    call addto(incompfix,q_c,scale=-1.0/rho_w)
    call addto(incompfix,q_r,scale=-1.0/rho_w)
    call scale(incompfix,density_local)
    call addto(incompfix,1.0)
        
    incompressible = ((gstat/=0).or.(cstat/=0))
    if(incompressible) then
       ewrite(0,*) "Selected compressible eos but not specified either C_P or C_V."
    end if
        
    call zero(drhodp)
    
    if(.not.incompressible) then
       energy_local=>extract_scalar_field(state,'InternalEnergy',stat=stat)
       
       if((gstat==0).and.(cstat==0)) then
          
          call allocate(energy_remap, drhodp%mesh, 'RemappedInternalEnergy')
          call remap_field(energy_local, energy_remap)
          
          do node=1,node_count(drhodp)
             ! Evil node computation goes here
             temperature_node=node_val(energy_remap,node)&
                  /node_val(c_v_bulk,node)
             drhodp_node=node_val(R_bulk,node)&
                  *temperature_node
             
             call set(drhodp, node, drhodp_node)
             if (present(temperature)) &
                  & call set(temp_local, node, temperature_node)
             ! horrible temporary hack
          end do

          if (has_scalar_field(state,"PotentialTemperature")) &
               call scale(drhodp,c_p) 
          call invert(drhodp)
          call scale(drhodp,incompfix)
          call scale(drhodp,incompfix)
          if (has_scalar_field(state,"PotentialTemperature")) &
               call scale(drhodp,c_v) 

       end if
       
       
    end if

    if(present(density)) then
       ! calculate the density
       ! density may equal density in state depending on how this
       ! subroutine is called
       if(incompressible) then
          ! density = reference_density
          call set(density, reference_density)
       else

          call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
          call remap_field(complete_pressure, pressure_remap)
!          call set(density, pressure_remap)
!          call scale(density, drhodp)


          do node=1,node_count(drhodp)
             ! Evil node computation goes here
             pressure_node=node_val(pressure_remap,node)
             temperature_node=node_val(energy_remap,node)&
                  /node_val(c_v_bulk,node)
             R_node=node_val(R_bulk,node)

             i_node=(node_val(q_r,node)+node_val(q_c,node))/rho_w

             density_node=pressure_node/&
                  (R_node*temperature_node+pressure_node*i_node)
             
             call set(density, node, density_node)
          end do

!          call set(density,energy_remap)
!          call scale(density,R_bulk)
!          call invert(density)
!          call scale(density,c_v_bulk)
!          call scale(density,pressure_remap)
          
!          call zero(incompfix)
!          call addto(incompfix,q_r,scale=1.0/rho_w)
!          call addto(incompfix,q_c,scale=1.0/rho_w)
!          call scale(incompfix,density)
!          call addto(incompfix,q_g)
!          call invert(incompfix)
!          call scale(density,incompfix)

!          call set_dirichlet_consistent(density)
          
          call deallocate(pressure_remap)
       end if
    end if
    
    if(present(pressure)) then
       if(incompressible) then
          ! pressure is unrelated to density in this case
          call zero(pressure)
       else
          ! calculate the pressure using the eos and the calculated (probably prognostic)
          ! density
          
          density_local=>extract_scalar_field(state,'Density',stat=stat)
          
          if (stat==0) then
             assert(pressure%mesh==drhodp%mesh)
             
             call allocate(density_remap, drhodp%mesh, "RemappedDensity")
             call remap_field(density_local, density_remap)


             do node=1,node_count(drhodp)
                ! Evil node computation goes here
                temperature_node=node_val(energy_remap,node)&
                     /node_val(c_v_bulk,node)
                R_node=node_val(R_bulk,node)
                i_node=(node_val(q_r,node)+node_val(q_c,node))/rho_w
                density_node=node_val(density_remap,node)

                pressure_node=density_node*R_node*temperature_node&
                     /(1.0-density_node*i_node)

                call set(pressure, node, pressure_node)
             end do

!             call set(pressure, drhodp)
!             call invert(pressure)
!             call scale(pressure, density_remap)
             
             call deallocate(density_remap)
          else
             FLExit('No Density in material_phase::'//trim(state%name))
          end if
       end if
    end if
    
    call deallocate(R_bulk)
    call deallocate(c_v_bulk)
    call deallocate(q_g)
    call deallocate(incompfix)
    if(.not.incompressible) call deallocate(energy_remap)
    call deallocate(dummyscalar)
    if (present(temperature)) then
       call set(temperature,temp_local)
       call deallocate(temp_local)
    end if
    
  end subroutine compressible_eos_giraldo_1mat
  
  subroutine compressible_eos_foam(state, eos_path, drhodp, &
    density, pressure)
    ! Foam EoS Used with compressible simulations of liquid drainage in foams.
    ! It describes the liquid content in the foam as the product of the  Plateau 
    ! border cross sectional area and the local Plateau  border length per unit volume (lambda).
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure

    ! locals
    integer :: pstat, dstat
    type(scalar_field), pointer :: pressure_local, density_local, drainagelambda_local
    real :: atmospheric_pressure
    type(scalar_field) :: pressure_remap, density_remap, drainagelambda_remap

    call zero(drhodp)

    pressure_local => extract_scalar_field(state,'Pressure', stat=pstat)

    drainagelambda_local => extract_scalar_field(state,'DrainageLambda')

    call allocate(drainagelambda_remap, drhodp%mesh, 'RemappedDrainageLambda')
    call remap_field(drainagelambda_local, drainagelambda_remap)

    call addto(drhodp, drainagelambda_remap)

    call deallocate(drainagelambda_remap)

    if(present(density)) then
      if (pstat==0) then
        assert(density%mesh==drhodp%mesh)

        call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
        call remap_field(pressure_local, pressure_remap)

        call set(density, atmospheric_pressure)
        call addto(density, pressure_remap)
        call scale(density, drhodp)

        call deallocate(pressure_remap)
      else
        FLExit('No Pressure in material_phase::'//trim(state%name))
      end if
    end if

    if(present(pressure)) then
      density_local=>extract_scalar_field(state,'Density',stat=dstat)
      if (dstat==0) then
        assert(pressure%mesh==drhodp%mesh)

        call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call allocate(density_remap, drhodp%mesh, "RemappedDensity")
        call remap_field(density_local, density_remap)

        call set(pressure, drhodp)
        call invert(pressure)
        call scale(pressure, density_remap)

        call deallocate(density_remap)
      else
        FLExit('No Density in material_phase::'//trim(state%name))
      end if
    end if
        

  end subroutine compressible_eos_foam
        
  subroutine compressible_material_eos(state,materialdensity,&
                                    materialpressure,materialdrhodp)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional :: materialdensity, &
                                             materialpressure, materialdrhodp

    !locals
    integer :: stat, gstat, cstat
    type(scalar_field), pointer :: pressure, materialenergy, materialdensity_local
    character(len=4000) :: thismaterial_phase, eos_path
    real :: reference_density, ratio_specific_heats
    real :: bulk_sound_speed_squared, atmospheric_pressure
    type(scalar_field) :: drhodp

    ewrite(1,*) 'Entering compressible_material_eos'

    if (present(materialdensity)) then
      call allocate(drhodp, materialdensity%mesh, 'Gradient of density wrt pressure')
    else if (present(materialpressure)) then
      call allocate(drhodp, materialpressure%mesh, 'Gradient of density wrt pressure')
    else if (present(materialdrhodp)) then
      call allocate(drhodp, materialdrhodp%mesh, 'Gradient of density wrt pressure')
    else
      FLAbort("No point in being in here if you don't want anything out.")
    end if

    thismaterial_phase = '/material_phase::'//trim(state%name)
    eos_path = trim(thismaterial_phase)//'/equation_of_state'

    if(have_option(trim(eos_path)//'/compressible')) then

      if(have_option(trim(eos_path)//'/compressible/stiffened_gas')) then
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/reference_density', &
                        reference_density, default=0.0)
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/ratio_specific_heats', &
                        ratio_specific_heats, stat=gstat)
        if(gstat/=0) then
          ratio_specific_heats=1.0
        end if
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/bulk_sound_speed_squared', &
                        bulk_sound_speed_squared, stat = cstat)
        if(cstat/=0) then
          bulk_sound_speed_squared=0.0
        end if
        if((gstat/=0).and.(cstat/=0)) then
          FLExit("Must set either a bulk_sound_speed_squared or a ratio_specific_heats.")
        end if
        materialenergy=>extract_scalar_field(state,'MaterialInternalEnergy',stat=stat)
        if(stat==0) then   ! we have an internal energy field
          drhodp%val = 1.0/( bulk_sound_speed_squared + (ratio_specific_heats - 1.0)*materialenergy%val )
        else               ! we don't have an internal energy field
          call set(drhodp, 1.0/bulk_sound_speed_squared)
        end if

        if(present(materialdensity)) then
          ! calculate the materialdensity
          ! materialdensity can equal materialdensity in state depending on how this
          ! subroutine is called
          pressure=>extract_scalar_field(state,'Pressure',stat=stat)
          if (stat==0) then
            call get_option(trim(pressure%option_path)//'/prognostic/atmospheric_pressure', &
                            atmospheric_pressure, default=0.0)
            materialdensity%val = drhodp%val*(pressure%val + atmospheric_pressure &
                                   + bulk_sound_speed_squared*reference_density)
          else
            FLExit('No Pressure in material_phase::'//trim(state%name))
          end if
        end if

        if(present(materialpressure)) then
          ! calculate the materialpressure using the eos and the calculated (probably prognostic)
          ! materialdensity
          ! materialpressure /= bulk pressure
          materialdensity_local=>extract_scalar_field(state,'MaterialDensity',stat=stat)
          if (stat==0) then
            materialpressure%val = materialdensity_local%val/drhodp%val &
                                  - bulk_sound_speed_squared*reference_density
          else
            FLExit('No MaterialDensity in material_phase::'//trim(state%name))
          end if
        end if

!       else
!       ! place other compressible material eos here

      end if

!     else
!     ! an incompressible option?

    end if      

    if(present(materialdensity)) then
      ewrite_minmax(materialdensity)
    end if

    if(present(materialpressure)) then
      ewrite_minmax(materialpressure)
    end if

    if(present(materialdrhodp)) then
      materialdrhodp%val=drhodp%val
      ewrite_minmax(materialdrhodp)
    end if

    call deallocate(drhodp)

  end subroutine compressible_material_eos

  subroutine subtract_hydrostatic_pressure_contribution(state,pressure,hp)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: pressure
    type(scalar_field), intent(in) :: hp

    type(scalar_field) :: rhs

    integer :: ele

    if ((continuity(pressure)<= continuity(hp)) .and. &
        (element_degree(pressure,1) >= element_degree(hp,1))) then

       ! No reason to make life difficult unnecessarily
       call addto(pressure,hp,scale=-1.0)
    else
       
       call allocate(rhs,pressure%mesh,"PerturbationPressureRHS")

       element_loop: do ele=1,ele_count(pressure)
          call assemble_perturbation_pressure_rhs(ele,pressure,hp,rhs)
       end do element_loop

!       call petsc_solve(pressure,get_mass_matrix(state,pressure%mesh),rhs,state)

       call deallocate(rhs)
    end if

  end subroutine subtract_hydrostatic_pressure_contribution
          
  subroutine assemble_perturbation_pressure_rhs(ele,pressure,hp,rhs)
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: pressure,hp
    type(scalar_field),intent(inout) :: rhs

    

  end subroutine assemble_perturbation_pressure_rhs


  subroutine set_temperature_from_potential_temperature(theta,pressure,&
       temperature,p0,R,c_p)
    type(scalar_field), intent(in) :: theta,pressure
    real, intent(in) :: p0,R,c_p
    type(scalar_field), intent(inout) :: temperature
    type(scalar_field) :: remap_theta, remap_pressure
    integer :: node

    call allocate(remap_theta,temperature%mesh,"RemappedTheta")
    call allocate(remap_pressure,temperature%mesh,"RemappedPressure")

    call remap_field(theta,remap_theta)
    call remap_field(pressure,remap_pressure)

    do node=1,node_count(temperature)

       call set(temperature,node,node_val(remap_theta,node)*(node_val(pressure,node)/p0)**(R/c_p))

    end do

    call deallocate(remap_pressure)
    call deallocate(remap_theta)

  end subroutine set_temperature_from_potential_temperature

    

  subroutine set_EOS_pressure_and_temperature(state)
    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), pointer :: temperature,nl_field, pressure, density
    integer :: work

    work=0

    ! Here there will be magic done to put the Equation of state pressure
    ! and insitu bulk Temperature on the density mesh in the input state

    if (has_scalar_field(state(1),"EOSPressure")) then
       pressure=>extract_scalar_field(state(1),"EOSPressure")
       work=work+1
    else
!       ensure
!       FLAbort("No EOSPressure field in state!")
    end if
    if (has_scalar_field(state(1),"InsituTemperature")) then
       temperature=>extract_scalar_field(state(1),"InsituTemperature")
       work=work+2
    else
!       FLAbort("No Temperature field in state!")
    end if

    select case(work)
    case(0)
    case(1)
       call compressible_eos(state,full_pressure=pressure)
    case(2)   
       call compressible_eos(state,temperature=temperature)
    case(3)   
       call compressible_eos(state,full_pressure=pressure,&
            temperature=temperature)
    end select
  contains 


      subroutine ensure_field(state,name)
        type(state_type) :: state
        character(len=*) :: name
        type(scalar_field) :: sfield1
        type(scalar_field), pointer :: sfield2
    
        if (.not. has_scalar_field(state,name)) then
           sfield2=>extract_scalar_field(state,"Density")
           call allocate(sfield1,sfield2%mesh,trim(name))
           call insert(state,sfield1,trim(sfield1%name))
           call zero(sfield1)
        endif

      end subroutine ensure_field

  end subroutine set_EOS_pressure_and_temperature

end module equation_of_state

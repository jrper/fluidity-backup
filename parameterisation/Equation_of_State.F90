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
  use sparsity_patterns_meshes
  use transform_elements, only: transform_to_physical
  use diagnostic_fields, only: calculate_galerkin_projection
  use vtk_interfaces, only: vtk_write_fields
  
  implicit none
  
  interface compressible_eos
     module procedure compressible_eos_1mat,&
          compressible_eos_mmat
  end interface compressible_eos

  private
  public :: calculate_perturbation_density, mcD_J_W_F2002, &
            compressible_eos, compressible_material_eos, &
            set_EOS_pressure_and_temperature, safe_set, &
            extract_entropy_variable

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
        
 
         complete_pressure=>extract_scalar_field(state,"Pressure")
 
         if (present(pressure))then
            call compressible_eos_giraldo_1mat(state, eos_path, drhodp_local, &
                 complete_pressure,density,pressure=pressure,temperature=temperature)
         else
            call compressible_eos_giraldo_1mat(state, eos_path, drhodp_local, &
                 complete_pressure,density=density, pressure=full_pressure,temperature=temperature)
         end if


      elseif(have_option(trim(eos_path)//'/compressible/foam')) then
        
        ! eos used in foam modelling
        
         call compressible_eos_foam(state, eos_path, drhodp_local, &
              density=density, pressure=pressure)

     end if
    
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

  subroutine compressible_eos_mmat(state, density, pressure,&
       drhodp,temperature,getold)

    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), intent(inout), optional :: density, drhodp
    type(scalar_field), intent(inout), optional :: pressure, temperature
    logical, intent(in), optional :: getold

    integer :: stat
    
    character(len=OPTION_PATH_LEN) :: eos_path
    logical :: getoldlocal = .false.

    if (size(state)==1) then
       call compressible_eos_1mat(state(1), density, pressure,&
            drhodp,temperature)
       return
    end if
    
    if (present(getold)) then
       getoldlocal=getold
    else
       getoldlocal=.false.
    end if

    ewrite(1,*) 'Entering compressible_eos_mmat'
    
    if (present(drhodp)) then
       if (present(density)) then
          assert(drhodp%mesh==density%mesh)
       end if
       if (present(pressure)) then
          assert(drhodp%mesh==pressure%mesh)
       end if
    end if
    if (.not. (present(drhodp) .or. present(density) .or.&
         present(pressure) .or. present(temperature))) then
       FLAbort("No point in being in here if you don't want anything out.")
    end if
   
    eos_path = trim(state(1)%option_path)//'/equation_of_state'
    
    if(have_option(trim(eos_path)//'/compressible')) then       
       if(have_option(trim(eos_path)//'/compressible/ATHAM')) then
          
          !ATHAM equation of state
          if (present(pressure)) then
             if (present(drhodp)) then
                call compressible_eos_atham(state,pressure,"Pressure",&
                     output2=drhodp,getold=getoldlocal)
             else
                call compressible_eos_atham(state,pressure,&
                     "Pressure",getold=getoldlocal)
             end if
          else if (present(density)) then
             if (present(drhodp)) then
                call compressible_eos_atham(state,density,"Density",&
                     output2=drhodp,getold=getoldlocal)
             else
                call compressible_eos_atham(state,density,&
                     "Density",getold=getoldlocal)
             end if
          else if (present(temperature)) then
             call compressible_eos_atham(state,temperature,&
                  "Temperature",getold=getoldlocal)
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

    if(present(drhodp)) then      
      ewrite_minmax(drhodp)
   end if

   if(present(temperature)) then      
      ewrite_minmax(temperature)
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
    
  end subroutine make_bulk_quantities_mmat

  subroutine safe_set(state,field1,field2)
    type(state_type) , intent(inout) :: state
    type(scalar_field), intent(inout) :: field1
    type(scalar_field), intent(in) :: field2

    if (continuity(field1)<=continuity(field2) .and. &
       element_degree(field1,1)>=element_degree(field2,1)) then
       call zero(field1)
       call remap_field(field2,field1)
    else
!       call zero(field1)
       call calculate_galerkin_projection(state,field2,field1,alpha=0.0)
       ewrite_minmax(field1)     
       ewrite_minmax(field2)
       call vtk_write_fields(field1%name,&
            position=extract_vector_field(state,&
            "Coordinate"),&
            model=field2%mesh,sfields=(/field1,field2,&
            extract_scalar_field(state,"Time")/))
    end if




  end subroutine safe_set

    

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

  subroutine make_atham_quantities(state,qc_p,qc_v,qc_solid,scaledq,old)

    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), intent(inout) :: qc_p,qc_v,qc_solid,scaledq
    character(len=*), intent(in) :: old

    type(scalar_field), pointer :: fraction, thermal
    type(scalar_field) :: qc_p_local,qc_v_local,&
         qc_solid_local,scaledq_local
    character(len=OPTION_PATH_LEN) :: eos_path
    integer :: i,stat
    real :: c_v,c_p,rho
    type(scalar_field) :: dry_fraction

    thermal=>extract_scalar_field(state(1),"InternalEnergy",stat=stat)
    if (stat /= 0) &
         thermal=>extract_scalar_field(state(1),"PotentialTemperature",stat=stat)
    if (stat /= 0) &
         thermal=>extract_scalar_field(state(1),"ATHAMInternalEnergy",stat=stat)
    if (stat /= 0) &
         FLAbort("Can't find thermodynamic variable in state!")


    call allocate(qc_p_local,thermal%mesh,"LocalqC_p")
    call allocate(qc_v_local,thermal%mesh,"LocalqC_v")
    call allocate(qc_solid_local,thermal%mesh,"LocalqC_solid")
    call allocate(scaledq_local,thermal%mesh,"LocalScaledQ")
    call allocate(dry_fraction,thermal%mesh,"DryGasMassFraction")


    call zero(qc_p_local)
    call zero(qc_v_local)
    call zero(qc_solid_local)
    call zero(scaledq_local)

    call set(dry_fraction,1.0)

    do i=2,size(state)
       if (has_scalar_field(state(i),"MassFraction")) then
          fraction=>extract_scalar_field(state(i),trim(old)//"MassFraction")
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
             call addto(qc_solid_local,fraction,scale=c_v)
             call addto(scaledq_local,fraction,scale=1.0/rho)
          else
             call addto(qc_v_local,fraction,scale=c_v)
             call addto(qc_p_local,fraction,scale=c_p)
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
    call addto(qc_v_local,dry_fraction,scale=c_v)
    call addto(qc_p_local,dry_fraction,scale=c_p)


    call set(qc_p,c_p)
    call set(qc_v,c_v)
    call set(qc_solid,0.0)
    call set(scaledq,0.0)

    call safe_set(state(1),qc_p,qc_p_local)
    call safe_set(state(1),qc_v,qc_v_local)
    call safe_set(state(1),qc_solid,qc_solid_local)
    call safe_set(state(1),scaledq,scaledq_local)

    call deallocate(dry_fraction)
    call deallocate(qc_p_local)
    call deallocate(qc_v_local)
    call deallocate(qc_solid_local)
    call deallocate(scaledq_local)
    
  end subroutine make_atham_quantities


  subroutine compressible_eos_atham(state,output,output_name,output2,getold)
    ! Routines to invert the ATHAM equation of state in various
    ! permuations
    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), intent(inout) :: output
    character(len=*), intent(in) :: output_name
    type(scalar_field), intent(inout), optional :: output2
    logical, intent(in), optional :: getold
    
    ! locals
    integer :: stat
    type(scalar_field), pointer :: pressure, thermal,density
    type(scalar_field) :: thermal_remap, pressure_remap, density_remap, &
                          & qc_v,qc_p,qc_solid,scaledq

    real :: p0, TV_node, rho_node,p_node,qc_p_node, &
         qc_v_node,qc_solid_node, &
         scaledq_node, drhodp_node,temp_node
    integer :: node, thermal_variable
    character (len=10) :: old="          "

    call get_option(trim(state(1)%option_path)//'/compressible/ATHAM/reference_pressure', &
                    p0, default=1.0e5)

    call allocate(qc_p,output%mesh,"qc_p")
    call allocate(qc_v,output%mesh,"qc_v")
    call allocate(qc_solid,output%mesh,"qc_solid")
    call allocate(scaledq,output%mesh,"scaledp")

    if (present(getold)) then
       if (getold) then
          old="Old"
       else
          old="   "
       end if
    else
       old="   "
    end if

    call make_atham_quantities(state,qc_p,qc_v,qc_solid,scaledq,trim(old))




    density=>extract_scalar_field(state(1),trim(old)//"Density",stat=stat)
    if (stat /=0 ) then
       FLAbort("No density in bulk state for ATHAM EoS")
    end if
    pressure=>extract_scalar_field(state(1),trim(old)//"Pressure",stat=stat)
    if (stat /=0 ) then
       FLAbort("No pressure in bulk state for ATHAM EoS")
    end if
    if (has_scalar_field(state(1),'InternalEnergy')) then
       thermal_variable=1
       thermal=>extract_scalar_field(state(1),trim(old)//'InternalEnergy',stat=stat)
    else if (has_scalar_field(state(1),'PotentialTemperature')) then
       thermal_variable=2
       thermal=>extract_scalar_field(state(1),trim(old)//'PotentialTemperature',&
            stat=stat)
    else if (has_scalar_field(state(1),'ATHAMInternalEnergy')) then
       thermal_variable=3
       thermal=>extract_scalar_field(state(1),trim(old)//'ATHAMInternalEnergy',&
            stat=stat)
    else
       FLAbort("No thermodynamic variable in bulk state for ATHAM EoS!")
    end if


    call allocate(thermal_remap,output%mesh,"Remeshed"//trim(thermal%name))
    call allocate(density_remap,output%mesh,"RemeshedDensity")
    call allocate(pressure_remap,output%mesh,"RemeshedPressure")

    call set(pressure_remap,p0)
    call set(thermal_remap,275.0)
    call set(density_remap,1.0)

    call safe_set(state(1),density_remap,density)
    call safe_set(state(1),pressure_remap,pressure)
    call safe_set(state(1),thermal_remap,thermal)

    if (output_name=="Pressure") &
         call pressure_atham(output,density_remap,thermal_remap)
    if (output_name=="Density") &
         call density_atham(output,pressure_remap,thermal_remap)
    if (output_name=="Temperature") &
         call temperature_atham(output,thermal_remap,pressure_remap)
    if (present(output2)) then
       if (output_name=="Pressure") then
          call drhodp_atham(output2,density_remap,pressure_remap,&
               thermal_remap)
       else if (output_name=="Density") then
          call drhodp_atham(output2,density_remap,pressure_remap,&
               thermal_remap)
       else
          call drhodp_atham(output2,density_remap,pressure_remap,&
               thermal_remap)
       end if
       ewrite_minmax(output2)
    end if

    call deallocate(qc_p)
    call deallocate(qc_v)
    call deallocate(qc_solid)
    call deallocate(scaledq)
    call deallocate(pressure_remap)
    call deallocate(density_remap)
    call deallocate(thermal_remap)

  contains

      subroutine pressure_atham(p,rho,TV)
        type(scalar_field), intent(inout) :: p
        type(scalar_field), intent(in):: rho,TV

        type(scalar_field) :: rhs

        real :: fp_node

        integer :: nlin

        select case(thermal_variable)
        case(1)
           do node=1,node_count(p)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)

              p_node=rho_node*(qc_p_node-qc_v_node)*TV_node&
                   /((1.0-rho_node*scaledq_node)*(qc_v_node+qc_solid_node))
              
              call set(p,node,p_node)
           end do
        case(2)
           
           call allocate(rhs,p%mesh,"RHS")

           do node=1,node_count(p)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)



              p_node=p0**(-(qc_p_node-qc_v_node)/qc_v_node)&
                      *(rho_node*(qc_p_node-qc_v_node)*TV_node)&
                      **(qc_p_node/qc_v_node)

              if (abs(scaledq_node)>1e-8) then


                 temp_node=rho_node*(qc_p_node-qc_v_node)&
                   *(qc_p_node+qc_solid_node)&
                   *TV_node/(1.0-rho_node*scaledq_node)
                 fp_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node)


                 do nlin=1,200
                    p_node=p_node&
                         -(fp_node-temp_node)&
                         /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                         *qc_v_node+qc_solid_node)
                    
                    fp_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                         *qc_p_node+qc_solid_node)
                    
                    if (sqrt((fp_node-temp_node)**2)<1.0e-10) then
                       exit
                    end if
                 
                 
                 end do
              end if
              call set(p,node,p_node)
              call set(rhs,node,temp_node)


           end do

!           call petsc_solve_nonlinear(p,func,dfuncdx)

           call deallocate(rhs)

        case(3)

           call allocate(rhs,p%mesh,"RHS")

           do node=1,node_count(p)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)



              p_node=p0**(-(qc_p_node-qc_v_node)/qc_v_node)&
                      *(rho_node*(qc_p_node-qc_v_node)*TV_node/(qc_p_node+qc_solid_node))&
                      **(qc_p_node/qc_v_node)

              if (abs(scaledq_node)>1e-8) then


                 temp_node=rho_node*(qc_p_node-qc_v_node)&
                   *TV_node/(1.0-rho_node*scaledq_node)
                 fp_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node)


                 do nlin=1,200
                    p_node=p_node&
                         -(fp_node-temp_node)&
                         /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                         *qc_v_node+qc_solid_node)
                    
                    fp_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                         *qc_p_node+qc_solid_node)
                    
                    if (sqrt((fp_node-temp_node)**2)<1.0e-10) then
                       exit
                    end if
                 
                 
                 end do
              end if
              call set(p,node,p_node)
              call set(rhs,node,temp_node)


           end do

!           call petsc_solve_nonlinear(p,func,dfuncdx)

           call deallocate(rhs)

        end select

      end subroutine pressure_atham
          
      subroutine func(x)
        type(scalar_field), intent(inout) :: x
            
        do node=1,node_count(x)
           
           TV_node=node_val(thermal_remap,node)
           rho_node=node_val(density_remap,node)
           qc_p_node=node_val(qc_p,node)
           qc_v_node=node_val(qc_v,node)
           qc_solid_node=node_val(qc_solid,node)
           scaledq_node=node_val(scaledq,node)

           p_node=node_val(x,node)
           
           p_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                *qc_p_node+qc_solid_node)-temp_node
           
           call set(x,node,p_node)
        end do
      end subroutine func

      subroutine dfuncdx(x,J)
        type(scalar_field) :: x
        type(csr_matrix) :: J
        type(csr_sparsity) :: sparsity
        
        if (.not. associated(J%val)) then
           call allocate(sparsity, node_count(x),&
                node_count(x),node_count(x), &
                diag=.true., name="JacobianSparsity")
           call allocate(J,sparsity)
           call deallocate(sparsity)
        end if
        
        do node=1,node_count(x)
           
           
           TV_node=node_val(thermal_remap,node)
           rho_node=node_val(density_remap,node)
           qc_p_node=node_val(qc_p,node)
           qc_v_node=node_val(qc_v,node)
           qc_solid_node=node_val(qc_solid,node)
           scaledq_node=node_val(scaledq,node)
           
           p_node=node_val(x,node)
           
           call set(J,node,node,((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                *qc_v_node+qc_solid_node))
        end do
        
      end subroutine dfuncdx

      subroutine density_atham(rho,p,TV)
        type(scalar_field), intent(inout) :: rho
        type(scalar_field), intent(in):: p,TV

        select case(thermal_variable)
        case(1)
           do node=1,node_count(rho)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)

              rho_node=(qc_v_node+qc_solid_node)*p_node&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +(qc_v_node+qc_solid_node)*p_node*scaledq_node)
              
              call set(rho,node,rho_node)
           end do

        case(2)
           do node=1,node_count(rho)
              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)

              rho_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node)&
                   /((qc_p_node-qc_v_node)*(qc_p_node+qc_solid_node)*TV_node&
                   +p_node*scaledq_node&
                   *((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node))
              
              call set(rho,node,rho_node)

           end do
        case(3)
           do node=1,node_count(rho)
              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)

              rho_node=p_node*((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node)&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +p_node*scaledq_node&
                   *((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)&
                      *qc_p_node+qc_solid_node))
              
              call set(rho,node,rho_node)

           end do

        end select

      end subroutine density_atham

      subroutine temperature_atham(T,TV,p)
        type(scalar_field), intent(inout) :: T
        type(scalar_field), intent(in):: p,TV

        select case(thermal_variable)
        case(1)

           call set(T,qc_v)
           call addto(T,qc_solid)
           call invert(T)
           call scale(T,TV)

        case(2)

           do node=1,node_count(T)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)

              temp_node=TV_node&
                   *(qc_p_node+qc_solid_node)&
                   /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)*qc_p_node&
                   +qc_solid_node)
                   

              call set(T,node,temp_node)
              
           end do
        case(3)

           do node=1,node_count(T)

              TV_node=node_val(TV,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)

              temp_node=TV_node&
                   /((p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)*qc_p_node&
                   +qc_solid_node)
                   

              call set(T,node,temp_node)
              
           end do

        end select

      end subroutine temperature_atham

      subroutine drhodp_atham(drhodp,rho,p,TV)
        type(scalar_field), intent(inout) :: drhodp
        type(scalar_field), intent(in):: rho,p,TV
              
         
        select case(thermal_variable)
        case(1)
           do node=1,node_count(drhodp)

              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)


              drhodp_node=(qc_v_node+qc_solid_node)&
                   *(qc_p_node-qc_v_node)*TV_node&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +(qc_v_node+qc_solid_node)*p_node*scaledq_node)**2

              call set(drhodp, node, drhodp_node)
           end do
        case(2)
           do node=1,node_count(drhodp)


              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)

              temp_node=(p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)
              drhodp_node=(temp_node*qc_v_node+qc_solid_node)*(qc_p_node-qc_v_node)&
                   *(qc_p_node+qc_solid_node)*TV_node&
                   /((qc_p_node-qc_v_node)*(qc_p_node+qc_solid_node)*TV_node&
                   +p_node*(temp_node*qc_p_node+qc_solid_node)*scaledq_node)**2

              call set(drhodp, node, drhodp_node)
           end do

        case(3)
           do node=1,node_count(drhodp)


              TV_node=node_val(TV,node)
              rho_node=node_val(rho,node)
              p_node=node_val(p,node)
              qc_p_node=node_val(qc_p,node)
              qc_v_node=node_val(qc_v,node)
              qc_solid_node=node_val(qc_solid,node)
              scaledq_node=node_val(scaledq,node)

              temp_node=(p0/p_node)**((qc_p_node-qc_v_node)/qc_p_node)
              drhodp_node=(temp_node*qc_v_node+qc_solid_node)*(qc_p_node-qc_v_node)&
                   *TV_node&
                   /((qc_p_node-qc_v_node)*TV_node&
                   +p_node*(temp_node*qc_p_node+qc_solid_node)*scaledq_node)**2

              call set(drhodp, node, drhodp_node)
           end do

        end select

      end subroutine drhodp_atham


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


  subroutine set_EOS_pressure_and_temperature(state,&
    have_temperature,&
    have_EOSPressure,&
    have_EOSDensity,extras)
    type(state_type), intent(inout), dimension(:) :: state
    type(scalar_field), pointer :: temperature,thermal, pressure, density
    integer :: i, stat
    logical :: have_temperature,have_EOSPressure,have_EOSDensity
    type(scalar_field), dimension(:), intent(inout), target, optional :: extras

    ! Here there will be magic done to put the Equation of state pressure
    ! and insitu bulk Temperature on the density mesh in the input state

    i=1
    if (.not. (have_temperature .and. have_EOSPressure&
         .and. have_EOSDensity) ) then
       thermal=>extract_scalar_field(state(1),"InternalEnergy",stat=stat)
       if (stat /= 0) &
            thermal=>extract_scalar_field(state(1),"PotentialTemperature",stat=stat)
       if (stat /= 0) &
            thermal=>extract_scalar_field(state(1),"ATHAMInternalEnergy",stat=stat)
       if (stat /= 0) &
            FLAbort("Can't find thermodynamic variable in state!")
    end if

    if (have_temperature) then
       temperature=>extract_scalar_field(state(1),"InsituTemperature",stat=stat)
       if (stat == 0) call compressible_eos(state,temperature=temperature)
    else
       call allocate(extras(i),thermal%mesh,"IteratedInsituTemperature")
       temperature=>extras(i)
       call zero(temperature)
       call compressible_eos(state,temperature=temperature)
       call allocate(extras(i+1),thermal%mesh,"OldInsituTemperature")
       temperature=>extras(i+1)
       call zero(temperature)
       call compressible_eos(state,temperature=temperature,getold=.true.)
       i=i+2
       stat=0
    end if

    if (have_EOSPressure) then
       pressure=>extract_scalar_field(state(1),"EOSPressure",stat)
       if (stat == 0)     call compressible_eos(state,pressure=pressure)
    else
       call allocate(extras(i),thermal%mesh,"IteratedEOSPressure")
       pressure=>extras(i)
       call zero(pressure)
       call compressible_eos(state,pressure=pressure)
       call allocate(extras(i+1),thermal%mesh,"OldEOSPressure")
       pressure=>extras(i+1)
       call zero(pressure)
       call compressible_eos(state,pressure=pressure,getold=.true.)
       i=i+2
    end if
    if (have_EOSDensity) then
       density=>extract_scalar_field(state(1),"IteratedEOSDensity",stat)
       if (stat == 0) call compressible_eos(state,density=density)
    else
       call allocate(extras(i),thermal%mesh,"IteratedEOSDensity")
       density=>extras(i)
       call zero(density)
       call compressible_eos(state,density=density)
       call allocate(extras(i+1),thermal%mesh,"OldEOSDensity")
       density=>extras(i+1)
       call zero(density)
       call compressible_eos(state,density=density,getold=.true.)
       i=i+2
    end if

  end subroutine set_EOS_pressure_and_temperature

  function extract_entropy_variable(states,prefix,suffix,stat) result(sfield)
    type(state_type), dimension(:) :: states
    character(len=*), optional :: prefix, suffix
    type(scalar_field), pointer :: sfield
    integer, optional :: stat


    character(len=OPTION_PATH_LEN) ::lprefix,lsuffix
    integer :: i,j,lstat

    character(len=*), dimension(4), parameter ::&
         entropy_names= (/ "PotentialTemperature",&
                           "InternalEnergy      ",&
                           "Temperature         ",&
                           "ATHAMInternalEnergy "/)

    ! Routine searches (in order) for one of the entropy names above 
    ! in the states given to it, then returns the field
                           
    if (present(prefix)) then
       lprefix=prefix
    else
       lprefix=""
    end if

     if (present(suffix)) then
       lsuffix=suffix
    else
       lsuffix=""
    end if

    do j=1,size(entropy_names)
       sfield=>extract_scalar_field(states,&
            trim(lprefix)//trim(entropy_names(j))//trim(lsuffix),stat=lstat)
       if (present(stat)) stat=lstat
       if (lstat==0) return
    end do
    
    if (.not. present(stat)) then 
       FLAbort('Failed to find entropy field when using Atham Equation of state')
    end if

  end function extract_entropy_variable

end module equation_of_state

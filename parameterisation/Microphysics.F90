!    Copyright (C) 2008 Imperial College London and others.
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

module microphysics
  !!< This module implements a warm cloud microphysics model in Fluidity

  use spud
  use state_module
  use fields
  use sparse_tools
  use fefields
  use fetools
  use boundary_conditions
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN,&
       PYTHON_FUNC_LEN, CURRENT_DEBUG_LEVEL
  use solvers
  use python_state
  use sparsity_patterns_meshes
  use equation_of_state

  implicit none

  type logic_array
     logical :: have_temp=.false.,&
         have_EOSDensity=.false.,&
         have_EOSPressure=.false.
  end type logic_array

  private
  public calculate_microphysics_forcings, initialise_microphysics,&
       store_microphysics_source, calculate_gas_density, &
       calculate_incompressible_pressure_correction

contains

  subroutine store_microphysics_source(state,name)
    type(state_type),intent(inout) :: state
    character(len=*) :: name
    type(scalar_field), pointer :: sfield1, sfield2
    integer :: stat

    sfield1=>extract_scalar_field(state,name,stat=stat)
    if (stat /= 0) then
       ewrite(3,*) ('No field '//name//' in state to store!')
       return
    end if
    
    sfield2=>extract_scalar_field(state,'Old'//name,stat=stat)
    if (stat /= 0) then
       ewrite(3,*) ('No old field '//name//' in state to store in!')
       return
    end if

    call set(sfield2,sfield1)

  end subroutine store_microphysics_source

  subroutine calculate_microphysics_forcings(state,current_time,dt,adapt)
    type(state_type), intent(inout) :: state(:)
    type(scalar_field), pointer :: sfield1,sfield2
    real, intent(in) :: current_time
    real, intent(in) :: dt
    logical, intent(in), optional :: adapt

    character(len=OPTION_PATH_LEN) :: prefix
    integer :: i
    logical :: microphysics_on, python_microphysics
    type(scalar_field), dimension (:), allocatable :: extras
    type(logic_array) :: logic
    
    if (have_option("/cloud_microphysics/microphysical_model")) then
       prefix="/cloud_microphysics/microphysical_model"
    end if
       
    microphysics_on=have_option("/cloud_microphysics")
    python_microphysics=.not.have_option("/cloud_microphysics/fortran_microphysics")

    if (microphysics_on) then
       if (.not. present(adapt)) then
          if (size(state)>1) then
             do i=2,size(state)
                if (have_option(trim(state(i)%option_path)&
                     //'/equation_of_state/compressible/ATHAM')) &
                     call store_microphysics_source(state(i),&
                     "MassFractionSource")
             end do
          else
             call store_microphysics_source(state(1),"CloudWaterFractionSource")
             call store_microphysics_source(state(1),"RainWaterFractionSource")
             call store_microphysics_source(state(1),"WaterVapourFractionSource")
          end if
          i=3
          if (has_scalar_field(state(1),"InsituTemperature")) then
             call store_microphysics_source(state(1),"InsituTemperature")
             logic%have_temp=.true.
             i=i-1
          end if
          if (has_scalar_field(state(1),"EOSPressure")) then
             call store_microphysics_source(state(1),"EOSPressure")
             logic%have_EOSPressure=.true.
             i=i-1
          end if
          if (has_scalar_field(state(1),"EOSDensity")) then
             call store_microphysics_source(state(1),"EOSDensity")
             logic%have_EOSDensity=.true.
             i=i-1
          end if
          allocate(extras(2*i))
       end if

       call set_EOS_pressure_and_temperature(state,&
            logic%have_temp,&
            logic%have_EOSPressure,&
            logic%have_EOSDensity,extras=extras)
       if (python_microphysics) then
          call calculate_microphysics_from_python(state,&
               prefix,current_time,dt,extras,logic)
       else
          call calculate_microphysics_from_fortran(state,&
               current_time,dt,extras)
       end if
       if (present(adapt)) then
          if (size(state)>1) then
             do i=2,size(state)
                call store_microphysics_source(state(i),"MassFractionSource")
             end do
          else
             call store_microphysics_source(state(1),"CloudWaterFractionSource")
             call store_microphysics_source(state(1),"RainWaterFractionSource")
             call store_microphysics_source(state(1),"WaterVapourFractionSource")
          end if
          if (logic%have_temp) &
          call store_microphysics_source(state(1),"InsituTemperature")
          if (logic%have_EOSPressure) &
               call store_microphysics_source(state(1),"EOSPressure")
          if (logic%have_EOSDensity) &
               call store_microphysics_source(state(1),"EOSDensity")
       end if

       do i=1,size(extras)
          call deallocate(extras(i))
       end do
       deallocate(extras)
    end if

  end subroutine calculate_microphysics_forcings

  subroutine calculate_microphysics_from_python(state,prefix,current_time,dt,extras,logic)
    ! Set microphysical  source terms from python.
    type(state_type),intent(inout), dimension(:) :: state
    character(len=*), intent(in) :: prefix
    character(len = 30) :: buffer
    real, intent(in) :: current_time
    real, intent(in) :: dt
    type(logic_array), intent(in) :: logic
    type(scalar_field), dimension(:), intent(in) :: extras

    character(len=*), parameter :: par_prefix="/cloud_microphysics/parameters"

    character(len=OPTION_PATH_LEN) :: option_path
    character(len=PYTHON_FUNC_LEN) :: pycode

    real :: v

    integer :: i

#ifdef HAVE_NUMPY

    call python_reset()
    call python_add_states(state)

    do i=1,size(extras)
       call python_add_field(extras(i),state(1))
    end do

    write(buffer,*) current_time
    call python_run_string("time="//trim(buffer))
    write(buffer,*) dt
    call python_run_string("dt="//trim(buffer))
    write(buffer,*) CURRENT_DEBUG_LEVEL
    call python_run_string("import __builtin__")
    call python_run_string("__builtin__.current_debug_level="//trim(buffer))

    call python_run_string("parameters={}")
    call python_run_string("parameters['c_v']={}")
    call python_run_string("parameters['c_p']={}")
    call python_run_string("parameters['rho']={}")
    call python_run_string("parameters['radius']={}")

    ! put gravity in parameters
    call get_option('/physical_parameters/gravity/magnitude',v,default=0.0)
    call set_parameter_real("parameters","g",v)

    ! put reference pressure in parameters

    option_path=trim(state(1)%option_path)//'/equation_of_state/&
         &compressible/ATHAM/reference_pressure'
    call get_option(trim(option_path),v,default=1.0e5)
    call set_parameter_real("parameters","p0",v)

    ! other parameters

    call get_option(trim(par_prefix)//'/L_v',v,default=2.5e6)
    call set_parameter_real("parameters","L_v",v)



    do i=1,size(state)
       option_path=trim(state(i)%option_path)//'/equation_of_state/&
         &compressible/ATHAM/'
       call get_option(trim(option_path)//'C_V',v,default=0.0)
       call set_parameter_real("parameters['c_v']",state(i)%name,v)
       call get_option(trim(option_path)//'C_P',v,default=0.0)
       call set_parameter_real("parameters['c_p']",state(i)%name,v)
       if (have_option((trim(option_path)//'density'))) then
          call get_option(trim(option_path)//'density',v)
          call set_parameter_real("parameters['rho']",state(i)%name,v)
       end if
       if (have_option((trim(option_path)//'radius'))) then
          call get_option(trim(option_path)//'radius',v)
          call set_parameter_real("parameters['radius']",state(i)%name,v)
       end if
    end do
    
    
    if (have_option(trim(prefix))) then
       ! read the users code from flml
       call get_option(trim(prefix),pycode)
    else
       ! default option saves time, but requires default setup
       pycode="import fluidity.cloud_microphysics as cm"//NEW_LINE('A')//"cm.testing(states,dt,parameters)"
    end if
! And finally run the user's code
    call python_run_string(trim(pycode))    
#else
    ewrite(-1,*) "When configuring, make sure Numpy is found"
    FLExit("Python microphysics requires Numpy")
#endif

  contains

    subroutine set_parameter_real(pyname,name,v)
      character(len=*), intent(in) :: pyname, name
      real, intent(in) :: v
      character(len = 30) :: buffer

      write(buffer,*) v
      call python_run_string(trim(pyname)//'["'//trim(name)//'"]='//trim(buffer))

    end subroutine set_parameter_real

  end subroutine calculate_microphysics_from_python

subroutine calculate_incompressible_pressure_correction(state,pressure_correction)
  type(state_type), intent(in) :: state
  type(scalar_field), intent(inout) :: pressure_correction
  type(scalar_field), pointer :: density, q_c,q_r
  type(scalar_field) :: q_g
  real :: rho_w=1020.0

  density=>extract_scalar_field(state,"OldDensity")
  q_c=>extract_scalar_field(state,"OldCloudWaterFraction")
  q_r=>extract_scalar_field(state,"OldRainWaterFraction")

  call zero(pressure_correction)
  call addto(pressure_correction,q_c,scale=-1.0/rho_w)
  call addto(pressure_correction,q_r,scale=-1.0/rho_w)
  call scale(pressure_correction,density)
  call addto(pressure_correction,1.0)

end subroutine calculate_incompressible_pressure_correction

subroutine calculate_microphysics_from_fortran(state,current_time,dt,extras)
    ! Set microphysical  source terms from python.
    type(state_type),intent(inout), dimension(:) :: state
    real, intent(in) :: current_time
    real, intent(in) :: dt
    type(scalar_field), dimension(:), intent(in) :: extras

    integer :: i,ele

    type microphysics_field
       type(scalar_field), dimension(3) :: data
       type(scalar_field), pointer :: forcing, sinking_velocity
       logical :: has_forcing,has_sinking_velocity
    end type microphysics_field

    type(microphysics_field) :: pq_v,pq_r,pq_c,pp,prho,pS
    type(mesh_type), pointer :: mesh

    ! data arrays for the pass into the fortran routine
    ! These have the following order

    ! (1) Old Timelevel Value (input, projected)
    ! (2) Previous nonlinear iteration (input, projected)
    ! (3) Current nonlinear iteration (input, projected)
    ! (4) Sinking Velocity (output, on Microphysics mesh)
    ! (5) Microphysics forcing (output, on Microphysics mesh)


    real, dimension(5) :: q_v,q_r,q_c,p,rho,S
        

#ifdef FORTRAN_MICROPHYSICS
#ifndef FORTRAN_MICROPHYSICS_EXAMPLE
    interface 
       subroutine microphysics_main(time,timestep,lq_v,lq_r,lq_c,lp,lrho,lS)
         implicit none
         real :: time, timestep
         real, dimension(5) :: lq_v,lq_r,lq_c,lp,lrho,lS
       end subroutine microphysics_main
    end interface
#endif

    mesh=>extract_mesh(state(1),"MicrophysicsMesh")


    call extract_and_project(state,pq_v, "WaterVapour","MassFraction")
    call extract_and_project(state,pq_c, "CloudWater", "MassFraction")
    call extract_and_project(state,pq_r, "RainWater",  "MassFraction")
    call extract_and_project(state,pp,   "Bulk",       "Pressure")
    call extract_and_project(state,prho, "Bulk",       "Density")
    call extract_and_project(state,pS,   "Bulk",       "Entropy",&
         entropy=.true.)

    do ele=1,node_count(pp%data(1))

       ! Get element arrays

       call set_local_array(q_v,pq_v,ele)
       call set_local_array(q_c,pq_c,ele)
       call set_local_array(q_r,pq_r,ele)
       call set_local_array(p,pp,ele)
       call set_local_array(rho,prho,ele)
       call set_local_array(S,pS,ele)

       !   Call external routines

       call microphysics_main(current_time,dt,&
       q_v,q_r,q_c,p,rho,S)

       !   Call use results
       call store_result(pq_v,q_v,ele)
       call store_result(pq_c,q_c,ele)
       call store_result(pq_r,q_r,ele)
       call store_result(prho,rho,ele)
       call store_result(pS,S,ele)
    end do

    call clean_up(pq_v)
    call clean_up(pq_c)
    call clean_up(pq_r)
    call clean_up(pp)
    call clean_up(prho)
    call clean_up(pS)

#else
    FLAbort("Attempting to call external Fortran Microphysics without a linked external microphysics routine")
#endif

  contains 

    subroutine extract_and_project(lstate,mfield,material,fname,&
        entropy)
      type(state_type), dimension(:), intent(inout) :: lstate
      type(microphysics_field) :: mfield
      type(scalar_field), pointer :: sfield
      type(vector_field), pointer :: X
      character(len=*), intent(in) :: material, fname

      logical, intent(in), optional :: entropy
      logical :: lentropy
      
      integer :: i,stat
      type(mesh_type), pointer :: lmesh

      character(len=*), dimension(3), parameter::&
           old=(/ "Old     ",&
                  "Iterated",&
                  "        " /)
      

      if (present(entropy)) then
         lentropy=entropy
      else
         lentropy=.false.
      end if

      lmesh=>extract_mesh(lstate(1),"MicrophysicsMesh")

      X=>extract_vector_field(lstate(1),"Coordinate")

      do i=1,size(mfield%data)
         call allocate(mfield%data(i),lmesh,trim(old(i))//material//fname)
         if (lentropy) then
            sfield=>extract_entropy_variable(lstate,prefix=trim(old(i)))
         else
            sfield=>extract_scalar_field(extract_state(lstate,material),&
                 trim(old(i))//fname)
         end if
         call project_field(sfield,mfield%data(i),X)
      end do

      mfield%has_forcing=.false.
      mfield%has_sinking_velocity=.false.

      if (lentropy) then
         mfield%forcing=>extract_entropy_variable(state,&
              suffix="MicrophysicsSource",stat=stat)
         if (stat==0) mfield%has_forcing=.true.
         mfield%sinking_velocity=>extract_entropy_variable(state,&
              suffix="SinkingVelocity",stat=stat)
         if (stat==0) mfield%has_sinking_velocity=.true.
      else
         mfield%forcing=>extract_scalar_field(extract_state(state,material),&
              fname// "MicrophysicsSource",stat=stat)
         if (stat==0) mfield%has_forcing=.true.
         mfield%sinking_velocity=>extract_scalar_field(extract_state(state,material),&
                 fname//"SinkingVelocity",stat=stat)
            if (stat==0) mfield%has_sinking_velocity=.true.
      end if

    end subroutine extract_and_project

    subroutine clean_up(field)
      type(microphysics_field), intent(inout) :: field
      integer :: i

      do i=1,size(field%data)
         call deallocate(field%data(i))
      end do

      nullify(field%forcing)
      nullify(field%sinking_velocity)
         
    end subroutine clean_up


    subroutine set_local_array(local_array,fields,lele)
      real, dimension(:), intent(inout) :: local_array
      type(microphysics_field), intent(in) :: fields
      integer, intent(in) :: lele
      integer :: i
      
      local_array=0.0
      do i=1,size(fields%data)
         local_array(i)=node_val(fields%data(i),lele)
      end do
      
    end subroutine set_local_array

    subroutine store_result(field,vloc,n)
      type(microphysics_field), intent(inout) :: field
      real, intent(in), dimension(5) :: vloc
      integer, intent(in) :: n

      if (field%has_sinking_velocity)&
           call set(field%sinking_velocity,n,vloc(4))
      if (field%has_forcing)&
           call set(field%forcing,n,vloc(5))

    end subroutine store_result

end subroutine calculate_microphysics_from_fortran

#ifdef FORTRAN_MICROPHYSICS_EXAMPLE

subroutine microphysics_main(time,timestep,lq_v,lq_r,lq_c,lp,lrho,lS)
  implicit none
  real :: time, timestep
  real, dimension(5) :: lq_v,lq_r,lq_c,lp,lrho,lS
  
  real T,q_sat,dq
  real, parameter :: cv=1000.0, cp=714.,&
       c_w=4200.0,cv_v=1840.,&
       cp_v=1380., p0=10000.0,&
       L_v=2257000.0
  
  
  ! data arrays for the pass into the fortran routine
  ! These have the following order
  
  ! (1) Old Timelevel Value (input, projected)
  ! (2) Previous nonlinear iteration (input, projected)
  ! (3) Current nonlinear iteration (input, projected)
  ! (4) Sinking Velocity (output, on Microphysics mesh)
  ! (5) Microphysics forcing (output, on Microphysics mesh)
  
  T=lS(1)/(cp+(c_w-cp)*(lq_c(1)+lq_r(1))+(cp_v-cp)*lq_v(1))&
             *(lp(1)/p0)**((cp*(1.0-lq_c(1)-lq_r(1))+(cp_v-cp)*lq_v(1))&
             /(cv*(1.0-lq_c(1)-lq_r(1))+(cv_v-cv)*lq_v(1)))
  
  q_sat=lp(1)*0.18/(lrho(1)*(cp_v-cv_v)*T)
  
  if (lq_v(1)>q_sat) then
     dq=(lq_v(1)-q_sat)/timestep
  else
     dq=0
  end if
  
  lq_v(5)=-dq
  lq_c(5)=dq
  lq_r(5)=0.0

  lS(5) =L_v*dq

  


end subroutine microphysics_main

#endif


subroutine calculate_gas_density(state,gas_density)
  type(state_type), intent(in) :: state
  type(scalar_field), pointer :: density, q_c,q_r
  type(scalar_field) :: q_g, gas_density
  real :: rho_w=1020.0

  density=>extract_scalar_field(state,"Density")
  q_c=>extract_scalar_field(state,"CloudWaterFraction")
  q_r=>extract_scalar_field(state,"RainWaterFraction")

  !call allocate(q_g,q_c%mesh,"GasMassFraction")
  !call set(q_g,1.0)
  !call addto(q_g,q_r,scale=-1.0)
  !call addto(q_g,q_c,scale=-1.0)


  call zero(gas_density)
  call addto(gas_density,q_c,scale=-1.0/rho_w)
  call addto(gas_density,q_r,scale=-1.0/rho_w)
  call scale(gas_density,density)
  call addto(gas_density,1.0)
  call invert(gas_density)

  call scale(gas_density,density)
  !!call scale(gas_density,q_g)
  !call deallocate(q_g)

end subroutine calculate_gas_density

subroutine initialise_microphysics(state,current_time,dt)
  type(state_type), intent(inout) :: state(:)
  type(scalar_field), pointer :: gas_density
  real, intent(in) :: current_time
  real, intent(in) :: dt
  
  logical :: microphysics_on
  integer :: i

  ! In case it's needed later, so we can ensure everything we need exists
  
  microphysics_on=have_option('/cloud_microphysics')

  if (microphysics_on) then
     call set_EOS_pressure_and_temperature(state,.true.,.true.,.true.)
     call calculate_microphysics_forcings(state,current_time,dt)
  end if

contains

  subroutine ensure_field(state,name)
    type(state_type) :: state
    character(len=*) :: name
    type(scalar_field) :: sfield1
    type(scalar_field), pointer :: sfield2
   
!    sfield2=>extract_scalar_field(state,"Density")
!    call allocate(sfield1,sfield2%mesh,trim(name))
!    call insert(state,sfield1,trim(sfield1%name))
!    call zero(sfield1)

  end subroutine ensure_field

end subroutine initialise_microphysics


end module microphysics

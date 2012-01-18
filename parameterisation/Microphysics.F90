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
  use fetools
  use boundary_conditions
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN,&
       PYTHON_FUNC_LEN, CURRENT_DEBUG_LEVEL
  use solvers
  use python_state
  use sparsity_patterns_meshes
  use equation_of_state

  implicit none

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
       FLAbort('No field '//name//' in state to store!')
    end if
    
    sfield2=>extract_scalar_field(state,'Old'//name,stat=stat)
    if (stat /= 0) then
       FLAbort('No old field '//name//' in state to store in!')
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
    logical :: microphysics_on
    
    if (have_option("/cloud_microphysics/microphysical_model")) then
       prefix="/cloud_microphysics/microphysical_mdeol"
    end if
       
    if (size(state)==1) then
       microphysics_on=has_scalar_field(state(1),"CloudWaterFraction")
    else
       microphysics_on=.false.
       do i = 1, size(state)
          if (trim(state(i)%name)=="CloudWater") &
               microphysics_on=.true.
       end do
    end if

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
          call store_microphysics_source(state(1),"InsituTemperature")
          call store_microphysics_source(state(1),"EOSPressure")
       end if

       call set_EOS_pressure_and_temperature(state)    
       call calculate_microphysics_from_python(state,&
            prefix,current_time,dt)
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
          call store_microphysics_source(state(1),"InsituTemperature")
          call store_microphysics_source(state(1),"EOSPressure")
       end if

    end if

  end subroutine calculate_microphysics_forcings

  subroutine calculate_microphysics_from_python(state,prefix,current_time,dt)
    ! Set microphysical  source terms from python.
    type(state_type),intent(inout), dimension(:) :: state
    character(len=*), intent(in) :: prefix
    character(len = 30) :: buffer
    real, intent(in) :: current_time
    real, intent(in) :: dt

    character(len=*), parameter :: par_prefix="/cloud_microphysics/parameters"

    character(len=OPTION_PATH_LEN) :: option_path
    character(len=PYTHON_FUNC_LEN) :: pycode

    real :: v

    integer :: i

#ifdef HAVE_NUMPY
    call python_reset()
    call python_add_states(state)
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
    ewrite(-1,*) "When configuring, make sure NumPy is found"
    FLExit("Python microphysics requires NumPy")
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
  if (size(state)==1) then
     microphysics_on=has_scalar_field(state(1),"CloudWaterFraction")
  else
     microphysics_on=.false.
     do i = 1, size(state)
        if (trim(state(i)%name)=="CloudWater") &
             microphysics_on=.true.
     end do
  end if

  if (microphysics_on) then
     call set_EOS_pressure_and_temperature(state)
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

!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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

module detectors
  !!< Wrapper module for all things detector related

  use detector_data_types
  use detector_tools
  use detector_parallel

  !! Stuff we still need...
  use state_module
  use fields
  use field_options
  use spud
  use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN
  use diagnostic_variables
  use integer_hash_table_module

  implicit none

  private

  public :: initialise_detectors, read_detector_move_options

  interface detector_field
     module procedure detector_field_scalar, detector_field_vector
  end interface

  character(len=OPTION_PATH_LEN), parameter :: rk_gs_path="/explicit_runge_kutta_guided_search"

contains

  function detector_field_scalar(sfield)
    !!< Return whether the supplied field should be included in the .detector file
    logical :: detector_field_scalar
    type(scalar_field), target, intent(in) :: sfield
    
    if (sfield%option_path=="".or.aliased(sfield)) then
       detector_field_scalar=.false.
    else
       detector_field_scalar = have_option(&
            trim(complete_field_path(sfield%option_path)) // &
            "/detectors/include_in_detectors")
    end if

  end function detector_field_scalar

  function detector_field_vector(vfield)
    !!< Return whether the supplied field should be included in the .detector file
    logical :: detector_field_vector
    type(vector_field), target, intent(in) :: vfield

    if (vfield%option_path=="".or.aliased(vfield)) then
       detector_field_vector=.false.
    else
       detector_field_vector = have_option(&
            trim(complete_field_path(vfield%option_path)) // &
            "/detectors/include_in_detectors")
    end if

  end function detector_field_vector

  subroutine read_detector_move_options(detector_list, options_path)
    ! Subroutine to allocate the detector parameters, 
    ! including RK stages and update vector
    type(detector_linked_list), intent(inout) :: detector_list
    character(len=*), intent(in) :: options_path

    integer :: i, j, k
    real, allocatable, dimension(:) :: stage_weights
    integer, dimension(2) :: option_rank

    if (have_option(trim(options_path))) then

       call get_option(trim(options_path)//"/subcycles",detector_list%n_subcycles)
       call get_option(trim(options_path)//"/search_tolerance",detector_list%search_tolerance)

       ! Forward Euler options
       if (have_option(trim(options_path)//"/forward_euler_guided_search")) then
          detector_list%velocity_advection = .true.
          detector_list%n_stages = 1
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          detector_list%timestep_weights = 1.0
       end if

       ! Parameters for classical Runge-Kutta
       if (have_option(trim(options_path)//"/rk4_guided_search")) then
          detector_list%velocity_advection = .true.
          detector_list%n_stages = 4
          allocate(stage_weights(detector_list%n_stages*(detector_list%n_stages-1)/2))
          stage_weights = (/0.5, 0., 0.5, 0., 0., 1./)
          allocate(detector_list%stage_matrix(detector_list%n_stages,detector_list%n_stages))
          detector_list%stage_matrix = 0.
          k = 0
          do i = 1, detector_list%n_stages
             do j = 1, detector_list%n_stages
                if(i>j) then
                   k = k + 1
                   detector_list%stage_matrix(i,j) = stage_weights(k)
                end if
             end do
          end do
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          detector_list%timestep_weights = (/ 1./6., 1./3., 1./3., 1./6. /)
       end if

       ! Generic Runge-Kutta options
       if (have_option(trim(options_path)//trim(rk_gs_path))) then
          detector_list%velocity_advection = .true.
          call get_option(trim(options_path)//trim(rk_gs_path)//"/n_stages",detector_list%n_stages)

          ! Allocate and read stage_matrix from options
          allocate(stage_weights(detector_list%n_stages*(detector_list%n_stages-1)/2))
          option_rank = option_shape(trim(options_path)//trim(rk_gs_path)//"/stage_weights")
          if (option_rank(2).ne.-1) then
             FLExit('Stage Array wrong rank')
          end if
          if (option_rank(1).ne.size(stage_weights)) then
             ewrite(-1,*) 'size expected was', size(stage_weights)
             ewrite(-1,*) 'size actually was', option_rank(1)
             FLExit('Stage Array wrong size')
          end if
          call get_option(trim(options_path)//trim(rk_gs_path)//"/stage_weights",stage_weights)
          allocate(detector_list%stage_matrix(detector_list%n_stages,detector_list%n_stages))
          detector_list%stage_matrix = 0.
          k = 0
          do i = 1, detector_list%n_stages
             do j = 1, detector_list%n_stages
                if(i>j) then
                   k = k + 1
                   detector_list%stage_matrix(i,j) = stage_weights(k)
                end if
             end do
          end do

          ! Allocate and read timestep_weights from options
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          option_rank = option_shape(trim(options_path)//trim(rk_gs_path)//"/timestep_weights")
          if (option_rank(2).ne.-1) then
             FLExit('Timestep Array wrong rank')
          end if
          if (option_rank(1).ne.size(detector_list%timestep_weights)) then
             FLExit('Timestep Array wrong size')
          end if
          call get_option(trim(options_path)//trim(rk_gs_path)//"/timestep_weights",detector_list%timestep_weights)
       end if

       ! Boundary reflection
       if (have_option(trim(options_path)//trim("/reflect_on_boundary"))) then
          detector_list%reflect_on_boundary=.true.
       end if

       if (have_option(trim(options_path)//trim("/parametric_guided_search"))) then
          detector_list%tracking_method = GUIDED_SEARCH_TRACKING
       elseif (have_option(trim(options_path)//trim("/pure_guided_search"))) then
          detector_list%tracking_method = PURE_GS
       elseif (have_option(trim(options_path)//trim("/geometric_tracking"))) then
          detector_list%tracking_method = GEOMETRIC_TRACKING
          if (have_option(trim(options_path)//"/geometric_tracking/periodic_mesh")) then
             detector_list%track_periodic = .true.
             call get_option(trim(options_path)//"/geometric_tracking/periodic_mesh/mesh/name", detector_list%periodic_tracking_mesh)
          end if
       else
          if (check_any_lagrangian(detector_list)) then
             ewrite(-1,*) "Found lagrangian detectors, but no tracking options"
             FLExit('No lagrangian particle tracking method specified')
          end if
       end if

    else
       if (check_any_lagrangian(detector_list)) then
          ewrite(-1,*) "Found lagrangian detectors, but no timstepping options"
          FLExit('No lagrangian timestepping specified')
       end if
    end if

  end subroutine read_detector_move_options

end module detectors

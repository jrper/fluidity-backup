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

module merge_points

  use adjacency_lists
  use data_structures
  use fields
  use fldebug
  use linked_lists
  use node_owner_finder
  use pickers
  use sparse_tools

  implicit none
  
  private
  
  public :: allocate, deallocate, merge_points_filter, get_filter_id

  interface allocate
     module procedure allocate_merge_points_filter
  end interface allocate

  interface deallocate
     module procedure deallocate_merge_points_filter
  end interface deallocate

  integer, parameter, public :: default_merge_levels = 3
    
  real, parameter, public :: default_merge_tolerance = 1.0e2 * epsilon(0.0)

  type merge_points_filter
     integer :: levels, dim
     real, dimension(:,:), pointer :: points
     integer, dimension(:), allocatable :: offsets
     integer, dimension(:), allocatable :: indices
     real, dimension(:,:), allocatable :: bounding_box
  end type merge_points_filter

  contains

    subroutine allocate_merge_points_filter(filter,points,mask,level)

      !!! Allocate a quick & dirty octree structure to help find duplicate points

      type(merge_points_filter), intent(out) :: filter
      real, dimension(:,:), target, intent(in) :: points
      logical, dimension(size(points,2)), intent(in), optional :: mask
      integer, intent(in), optional :: level
      type(ilist), dimension(:), allocatable :: buckets

      integer :: dim, llevel, node

      dim=size(points,1)

      if (present(level)) then
         llevel=level
      else
         llevel=default_merge_levels
      end if

       filter%levels=llevel
       filter%dim=dim
       filter%points=>points

       allocate(buckets(2**(dim*llevel)))
       allocate(filter%bounding_box(dim,2))

       filter%bounding_box(:,1)=minval(points,dim=2)-default_merge_tolerance
       filter%bounding_box(:,2)=maxval(points,dim=2)+default_merge_tolerance


       !!! Place the live input points in the octree bins
       do node=1,size(points,2)
          if (present(mask)) then
             if (.not. mask(node)) cycle
          end if
          call place_point(filter,buckets,node)
       end do

       !!! Save memory by linearising the data back into plain arrays
       call collapse_buckets(filter,buckets)

       deallocate(buckets)

     end subroutine allocate_merge_points_filter

     subroutine place_point(filter,buckets,node)
       type(merge_points_filter), intent(inout) :: filter
       type(ilist), dimension(:) :: buckets
       integer, intent(in) :: node
       
       integer :: bucket(filter%dim), dim
       logical :: on_boundary(filter%dim)
       real :: r,delta
       

       do dim=1,filter%dim
          r=filter%points(dim,node)-filter%bounding_box(dim,1)
          delta=(filter%bounding_box(dim,2)-filter%bounding_box(dim,1))/2**filter%levels
          bucket(dim)=min(2**filter%levels,floor(r/delta)+1)
          on_boundary(dim)=abs(r-bucket(dim)*delta)<default_merge_tolerance&
               .and. bucket(dim)<2**filter%levels
       end do


       !!! Now put the point in the relevant bucket
       call put_in_buckets(bucket,on_boundary,1)


       contains
         
         recursive subroutine put_in_buckets(bucket_list,on_bound,id)

           integer, dimension(:) :: bucket_list
           logical, dimension(:) :: on_bound
           integer :: id
           
           !!! if the point is within the tolerence of the neighbouring bucket, stick it in both.

           if (size(bucket_list)==0) then
              call insert(buckets(id),node)
           else
              call put_in_buckets(bucket_list(2:),on_bound(2:),&
                   (id-1)*2**filter%levels+bucket_list(1))
              if (on_bound(1))&
                   call put_in_buckets(bucket_list(2:),on_bound(2:),&
                   (id-1)*2**filter%levels+bucket_list(1)+1)
           end if
          
         end subroutine put_in_buckets

       end subroutine place_point

     subroutine collapse_buckets(filter,buckets)
       type(merge_points_filter), intent(inout) :: filter
       type(ilist), dimension(:) :: buckets

       integer :: i, length

       length=0
       allocate(filter%offsets(size(buckets)+1))
       filter%offsets(1)=1
       do i=1,size(buckets)
          length=length+buckets(i)%length
          filter%offsets(i+1)=length+1
       end do
       allocate(filter%indices(length))
       
       do i=1,size(buckets)
          filter%indices(filter%offsets(i):filter%offsets(i+1)-1)=&
               list2vector(buckets(i))
          call deallocate(buckets(i))
       end do

     end subroutine collapse_buckets

     function get_filter_id(filter,point)
       type(merge_points_filter), intent(inout) :: filter
       real, dimension(filter%dim) :: point

       integer :: get_filter_id

       integer :: bucket, dim, i
       real :: r,delta

       !!! The utility routine. If point has a twin inside filter return
       !!! the index within the original filter point array, otherwise
       !!! return -1. If the point has muliple copies, only the first is returned.

       get_filter_id=-1

       !!! If point isn't in the bounding box, it's got to be new.
       if (any(point<filter%bounding_box(:,1)&
            .or. point>filter%bounding_box(:,2))) return

       bucket=1


       !!! Same hashing algorithm as before, but now we don't need to worry
       !!! about distance from the bucket boundaries because we already took
       !!! care of that.
       do dim=1,filter%dim
          bucket=(bucket-1)*2**filter%levels
          r=point(dim)-filter%bounding_box(dim,1)
          delta=(filter%bounding_box(dim,2)-filter%bounding_box(dim,1))/2**filter%levels
          bucket=bucket+min(2**filter%levels,floor(r/delta)+1)
       end do

       do i=filter%offsets(bucket),filter%offsets(bucket+1)-1
          !!! loop i
          if (any(abs(point-filter%points(:,filter%indices(i)))>&
               default_merge_tolerance)) then
             cycle
          else
             !!! success, of a sort
             get_filter_id=filter%indices(i)
             return
          end if
       end do

     end function get_filter_id
     
     subroutine deallocate_merge_points_filter(filter,deallocate_points)
       type(merge_points_filter), intent(inout) :: filter
       logical, optional :: deallocate_points
       
       !!! clean up

       if (present_and_true(deallocate_points)) then
          deallocate(filter%points)
       end if
       deallocate(filter%offsets)
       deallocate(filter%indices)
       deallocate(filter%bounding_box)
         
     end subroutine deallocate_merge_points_filter

end module merge_points

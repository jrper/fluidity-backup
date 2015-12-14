subroutine test_merge_points_filter

  use merge_points
  use unittest_tools
  implicit none

  type(merge_points_filter) :: filter
  real, dimension(2,5) :: points
  logical, dimension(5) :: mask
  
  logical :: fail

  points(:,1)= [0,0]
  points(:,2)= [1,0]
  points(:,3)= [2.5,1.4]
  points(:,4)= [3.0,2.0]
  points(:,5)= [7.3,1.2]

  mask=.true.
  mask(4)=.false.

  call allocate(filter,points,mask)

  fail = get_filter_id(filter,[0.,0.]) /= 1
  call report_test("[find zero]", fail, .false., "Should be 1")

  fail = get_filter_id(filter,[2.5,1.4]) /= 3
  call report_test("[find filter point]", fail, .false., "Should give 3")

  fail = get_filter_id(filter,[8.0,2.6]) /= -1
  call report_test("[find non filter point]", fail, .false., "Should give -1")
 
  fail = get_filter_id(filter,[3.0,2.0]) /= -1
  call report_test("[mask_works]", fail, .false., "Should be -1 !")

  fail = get_filter_id(filter,[1000.0,-50.0]) /= -1
  call report_test("[Outside_bounding_box]", fail, .false., "Should be -1 !")

  call deallocate(filter)
end subroutine test_merge_points_filter

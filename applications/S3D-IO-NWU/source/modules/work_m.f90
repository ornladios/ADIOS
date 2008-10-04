!=========================================================================================
  module work_m
!=========================================================================================
! module for work arrays
!
! really, work arrays shouldn't be needed.  Each array can be declared locally and
! most compilers should be smart enough to free up the memory space once the routine
! has exited...

  implicit none
!-----------------------------------------------------------------------------------------
! real arrays

  real, allocatable :: work1_1(:,:,:)     !generic work array of scalar size
  real, allocatable :: work1_2(:,:,:)     !generic work array of scalar size

!-----------------------------------------------------------------------------------------
  contains
!=========================================================================================
  subroutine allocate_work_arrays(flag)
!=========================================================================================
! allocate work arrays
!-----------------------------------------------------------------------------------------
  use param_m, only : nx, ny, nz, n_spec

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed

  integer flag
!-----------------------------------------------------------------------------------------
! work arrays

  if(flag.eq.1) then

    allocate(work1_1(nx,ny,nz));         work1_1=0.0
    allocate(work1_2(nx,ny,nz));         work1_2=0.0

  elseif(flag.eq.-1) then

    deallocate(work1_1)
    deallocate(work1_2)

  endif
!-----------------------------------------------------------------------------------------
  return
  end subroutine allocate_work_arrays
!-----------------------------------------------------------------------------------------
  end module work_m




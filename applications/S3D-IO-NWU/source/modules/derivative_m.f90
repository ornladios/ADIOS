!=========================================================================================
  module derivative_m
!=========================================================================================
! module for derivative variables
!-----------------------------------------------------------------------------------------
! Change Record
! 18-May-2005 Evatt Hawkes modified the derivative routines to communicate through a 
!             contiguous array, instead of strided MPI_datatypes. This was done to fix
!             an issue with the infiniband switch, 
!             and should be helpful on other platforms too.
! 24-NOV-2004 Evatt Hawkes - adding neg_f_x etc for single plane in-plane derivative ops
! 
!-----------------------------------------------------------------------------------------

  implicit none
!-----------------------------------------------------------------------------------------
! integers

  integer :: ibound=3
  integer :: i_symmetry=0 ! set as non-symmetric by default
                          ! must modify all calls to derivatives to make symmetry work
                          ! this is a big job and anyone who does it should speak with SDM

  integer :: isym_x=0     ! set as non-symmetric by default
  integer :: isym_y=0     ! set as non-symmetric by default
  integer :: isym_z=0     ! set as non-symmetric by default

! ghost cell arrays

#ifdef COARRAYCOMMUNICATION
  real, allocatable :: neg_f_x(:,:,:)[:], pos_f_x(:,:,:)[:]
  real, allocatable :: neg_f_y(:,:,:)[:], pos_f_y(:,:,:)[:]
  real, allocatable :: neg_f_z(:,:,:)[:], pos_f_z(:,:,:)[:]
#else
  real, allocatable :: neg_f_x(:,:,:), pos_f_x(:,:,:)
  real, allocatable :: neg_f_y(:,:,:), pos_f_y(:,:,:)
  real, allocatable :: neg_f_z(:,:,:), pos_f_z(:,:,:)
! 18-MAY-2005 Evatt Hawkes
! infiniband workarounds
  real, allocatable :: neg_fs_x(:,:,:), pos_fs_x(:,:,:)
  real, allocatable :: neg_fs_y(:,:,:), pos_fs_y(:,:,:)
#endif

! 24-NOV-2004 Evatt Hawkes
! ghost cell arrays for in-plane derivative operations

  real, target, allocatable :: neg_f_x_xy(:,:), pos_f_x_xy(:,:)
  real, target, allocatable :: neg_f_y_xy(:,:), pos_f_y_xy(:,:)

  real, target, allocatable :: neg_f_x_xz(:,:), pos_f_x_xz(:,:)
  real, target, allocatable :: neg_f_z_xz(:,:), pos_f_z_xz(:,:)

  real, target, allocatable :: neg_f_y_yz(:,:), pos_f_y_yz(:,:)
  real, target, allocatable :: neg_f_z_yz(:,:), pos_f_z_yz(:,:)

!-----------------------------------------------------------------------------------------
  contains
!=========================================================================================
  subroutine allocate_derivative_arrays(flag)
!=========================================================================================
! allocate derivative arrays
!-----------------------------------------------------------------------------------------
  use param_m, only : nx, ny, nz, iorder

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed

  integer flag
!-----------------------------------------------------------------------------------------
! derivative arrays

  if(flag.eq.1) then

#ifdef COARRAYCOMMUNICATION
      allocate(neg_f_x(iorder/2,ny,nz)[0:*])
      allocate(pos_f_x(iorder/2,ny,nz)[0:*])

      allocate(neg_f_y(nx,iorder/2,nz)[0:*])
      allocate(pos_f_y(nx,iorder/2,nz)[0:*])

      allocate(neg_f_z(nx,ny,iorder/2)[0:*])
      allocate(pos_f_z(nx,ny,iorder/2)[0:*])
#else
      allocate(neg_f_x(iorder/2,ny,nz))
      allocate(pos_f_x(iorder/2,ny,nz))

      allocate(neg_f_y(nx,iorder/2,nz))
      allocate(pos_f_y(nx,iorder/2,nz))

      allocate(neg_f_z(nx,ny,iorder/2))
      allocate(pos_f_z(nx,ny,iorder/2))

!     18-MAY-2005 Evatt Hawkes
!     infiniband work-arounds

      allocate(neg_fs_x(iorder/2,ny,nz))
      allocate(pos_fs_x(iorder/2,ny,nz))

      allocate(neg_fs_y(nx,iorder/2,nz))
      allocate(pos_fs_y(nx,iorder/2,nz))

#endif

!     24-NOV-2004 Evatt Hawkes
!     ghost cell arrays for in-plane derivative operations

      allocate(neg_f_x_xy(iorder/2,ny))
      allocate(pos_f_x_xy(iorder/2,ny))

      allocate(neg_f_y_xy(nx,iorder/2))
      allocate(pos_f_y_xy(nx,iorder/2))

      allocate(neg_f_x_xz(iorder/2,nz))
      allocate(pos_f_x_xz(iorder/2,nz))

      allocate(neg_f_z_xz(nx,iorder/2))
      allocate(pos_f_z_xz(nx,iorder/2))

      allocate(neg_f_y_yz(iorder/2,nz))
      allocate(pos_f_y_yz(iorder/2,nz))

      allocate(neg_f_z_yz(ny,iorder/2))
      allocate(pos_f_z_yz(ny,iorder/2))

  elseif(flag.eq.-1) then

      deallocate(neg_f_x)
      deallocate(pos_f_x)

      deallocate(neg_f_y)
      deallocate(pos_f_y)

      deallocate(neg_f_z)
      deallocate(pos_f_z)

#ifndef COARRAYCOMMUNICATION
!     18-MAY-2005 Evatt Hawkes
!     infiniband work-arounds

      deallocate(neg_fs_x)
      deallocate(pos_fs_x)

      deallocate(neg_fs_y)
      deallocate(pos_fs_y)
#endif

!     24-NOV-2004 Evatt Hawkes
!     ghost cell arrays for in-plane derivative operations

      deallocate(neg_f_x_xy)
      deallocate(pos_f_x_xy)

      deallocate(neg_f_y_xy)
      deallocate(pos_f_y_xy)

      deallocate(neg_f_x_xz)
      deallocate(pos_f_x_xz)

      deallocate(neg_f_z_xz)
      deallocate(pos_f_z_xz)

      deallocate(neg_f_y_yz)
      deallocate(pos_f_y_yz)

      deallocate(neg_f_z_yz)
      deallocate(pos_f_z_yz)

  endif
!-----------------------------------------------------------------------------------------
  return
  end subroutine allocate_derivative_arrays
!=========================================================================================
  subroutine initialize_derivative(io)
!=========================================================================================
! initializes derivative variables
!-----------------------------------------------------------------------------------------
  use topology_m

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io
!-----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    write(io,*) 'initializing derivative module...'
    write(io,*)
  endif
!-----------------------------------------------------------------------------------------
! allocate arrays

  call allocate_derivative_arrays(1)
!-----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    call write_header(io,'-')
  endif
!-----------------------------------------------------------------------------------------
  end subroutine initialize_derivative
!-----------------------------------------------------------------------------------------
  end module derivative_m


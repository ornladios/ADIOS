!=========================================================================================
  module topology_m
!=========================================================================================
! module for topology variables
!-----------------------------------------------------------------------------------------
! Change Record
!
! 24-NOV-2004: Evatt Hawkes - adding MPI types for single plane in-plane derivative ops
!
!-----------------------------------------------------------------------------------------


  implicit none
!-----------------------------------------------------------------------------------------
  include 'mpif.h'
!-----------------------------------------------------------------------------------------
! integers

  integer npes          !total number of processors
  integer myid          !rank of local processor
  integer xpes          !number of processors in x-direction
  integer xid           !id of local processor in x-direction
  integer ypes          !number of processors in y-direction
  integer yid           !id of local processor in y-direction
  integer zpes          !number of processors in z-direction
  integer zid           !id of local processor in z-direction

  integer term_status   !termination status (SDM)

  integer neighbor(6)
  integer xcomm, ycomm, zcomm, gcomm
  integer :: yz_comm      ! YZ plane communicator
  integer :: yz_pes, yz_id
  integer :: xz_comm      ! XZ plane communicator
  integer :: xz_pes, xz_id
  integer :: xy_comm      ! XY plane communicator
  integer :: xy_pes, xy_id

  integer xrfft_type,  yrfft_type,  zrfft_type
  integer xcfft_type,  ycfft_type,  zcfft_type
  integer xrows_type,  yrows_type,  zrows_type
  integer fxrows_type, fyrows_type, fzrows_type
  integer xcs_type,    ycs_type,    zcs_type

  integer ierr, status(MPI_Status_size)

! Addition 24-NOV-2004 Evatt Hawkes
  integer xy_xrows_type,  xy_yrows_type
  integer xz_xrows_type,  xz_zrows_type
  integer yz_yrows_type,  yz_zrows_type

! wkliao: make these 3 global, to be used in MPI I/O setup
  integer mypx, mypy, mypz

!-----------------------------------------------------------------------------------------
  contains
!========================================================================================
  subroutine initialize_topology(io,nx,ny,nz,npx,npy,npz,iorder,iforder)
!========================================================================================
! routine initializes some MPI stuff and the Cartesian MPI grid
!----------------------------------------------------------------------------------------
  implicit none
!----------------------------------------------------------------------------------------
! declarations passed in

  integer io
  integer nx, ny, nz, npx, npy, npz, iorder, iforder

! local declarations
  
  integer tmp
! integer mypx, mypy, mypz

  logical qflag
!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    write(io,*) 'initializing topology module...'
    write(io,*)
  endif
!----------------------------------------------------------------------------------------
! set term_status=0 for all processors

  term_status=0
!----------------------------------------------------------------------------------------
!	initialize MPI stuff

  term_status=0

  call MPI_Initialized(qflag, ierr)

  if(.not.qflag) then
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)

!   create communicator duplicate for global calls

    call MPI_Comm_dup(MPI_COMM_WORLD, gcomm, ierr)

  endif
!----------------------------------------------------------------------------------------
! check for npes compatibility

  if(npx*npy*npz.ne.npes) then
    if(myid.eq.0) then
      write(io,1) npx, npy, npz, npes
    endif
    call terminate_run(io,0)  !must be called by all processors
  endif

1 format(' npx*npy*npz is not equal to npes, npx = ',  &
            i5,' npy = ', i5, ' npz = ', i5, ' npes = ', i5)
!----------------------------------------------------------------------------------------
!	initialize Cartesian grid

  mypz = myid/(npx*npy)
  mypx = mod(myid-(mypz*npx*npy), npx)
  mypy = (myid-(mypz*npx*npy))/npx

  tmp = mypx-1 
  if(tmp.lt.0) then
    neighbor(1) = -1
  else 
    neighbor(1) = myid-1
  endif

  tmp = mypx+1 
  if (tmp.gt.npx-1) then
    neighbor(2) = -1
  else 
    neighbor(2) = myid+1
  endif

  tmp = mypy-1 
  if (tmp.lt.0) then
    neighbor(3) = -1
  else 
    neighbor(3) = myid-npx
  endif

  tmp = mypy+1 
  if (tmp.gt.npy-1) then
    neighbor(4) = -1
  else 
    neighbor(4) = myid+npx
  endif

  tmp = mypz-1 
  if (tmp.lt.0) then
    neighbor(5) = -1
  else 
    neighbor(5) = myid-(npx*npy)
  endif

  tmp = mypz+1 
  if (tmp.gt.npz-1) then
    neighbor(6) = -1
  else 
    neighbor(6) = myid+(npx*npy)
  endif
!----------------------------------------------------------------------------------------
!	create communicators for the x, y, and z directions

  call MPI_Comm_split(gcomm, mypy+1000*mypz, myid, xcomm,ierr)
  call MPI_Comm_split(gcomm, mypx+1000*mypz, myid, ycomm,ierr)
  call MPI_Comm_split(gcomm, mypx+1000*mypy, myid, zcomm,ierr)

  call MPI_Comm_rank(xcomm, xid, ierr)
  call MPI_Comm_size(xcomm, xpes,ierr)

  call MPI_Comm_rank(ycomm, yid, ierr)
  call MPI_Comm_size(ycomm, ypes,ierr)

  call MPI_Comm_rank(zcomm, zid, ierr)
  call MPI_Comm_size(zcomm, zpes,ierr)
!----------------------------------------------------------------------------------------
! Create MPI Comminicators for boundary planes.  This is used in the Boundary conditions

  call MPI_Comm_split(gcomm, xid, myid, yz_comm, ierr)   ! YZ plane communicator
  call MPI_Comm_split(gcomm, yid, myid, xz_comm, ierr)   ! XZ plane communicator
  call MPI_Comm_split(gcomm, zid, myid, xy_comm, ierr)   ! XY plane communicator

  call MPI_Comm_rank(yz_comm, yz_id, ierr)
  call MPI_Comm_size(yz_comm, yz_pes, ierr)

  call MPI_Comm_rank(xz_comm, xz_id, ierr)
  call MPI_Comm_size(xz_comm, xz_pes, ierr)

  call MPI_Comm_rank(xy_comm, xy_id, ierr)
  call MPI_Comm_size(xy_comm, xy_pes, ierr)

!----------------------------------------------------------------------------------------
! real FFT pencils

  call MPI_Type_vector(nx, 1,     1, MPI_REAL8, xrfft_type, ierr)
  call MPI_Type_vector(ny, 1,    nx, MPI_REAL8, yrfft_type, ierr)
  call MPI_Type_vector(nz, 1, nx*ny, MPI_REAL8, zrfft_type, ierr)

  call MPI_Type_commit(xrfft_type,ierr)
  call MPI_Type_commit(yrfft_type,ierr)
  call MPI_Type_commit(zrfft_type,ierr)
!-----------------------------------------------------------------------------------------
! complex FFT pencils

  call MPI_Type_vector(nx , 2,     1, MPI_REAL8, xcfft_type, ierr)
  call MPI_Type_vector(ny , 2,    nx, MPI_REAL8, ycfft_type, ierr)
  call MPI_Type_vector(nz , 2, nx*ny, MPI_REAL8, zcfft_type, ierr)

  call MPI_Type_commit(xcfft_type,ierr)
  call MPI_Type_commit(ycfft_type,ierr)
  call MPI_Type_commit(zcfft_type,ierr)
!-----------------------------------------------------------------------------------------
! derivative ghost planes

  call MPI_Type_vector(ny*nz,     iorder/2,    nx, MPI_REAL8, xrows_type, ierr)
  call MPI_Type_vector(   nz,    nx*(iorder/2), nx*ny, MPI_REAL8, yrows_type, ierr)
  call MPI_Type_vector(    1, nx*ny*(iorder/2),     1, MPI_REAL8, zrows_type, ierr)

  call MPI_Type_commit(xrows_type,ierr)
  call MPI_Type_commit(yrows_type,ierr)
  call MPI_Type_commit(zrows_type,ierr)
!-----------------------------------------------------------------------------------------
! derivative ghost planes for single plane in plane derivatives
! Addition 24-NOV-2004 Evatt Hawkes

  call MPI_Type_vector(ny*1 ,         iorder/2,    nx, MPI_REAL8,xy_xrows_type, ierr)
  call MPI_Type_vector(   1 ,    nx*(iorder/2), nx*ny, MPI_REAL8,xy_yrows_type, ierr)

  call MPI_Type_vector( 1*nz,         iorder/2,    nx, MPI_REAL8,xz_xrows_type, ierr)
  call MPI_Type_vector(    1,  nx*1*(iorder/2),     1, MPI_REAL8,xz_zrows_type, ierr)

  call MPI_Type_vector(   nz,     1*(iorder/2),  1*ny, MPI_REAL8,yz_yrows_type, ierr)
  call MPI_Type_vector(    1,  1*ny*(iorder/2),     1, MPI_REAL8,yz_zrows_type, ierr)

  call MPI_Type_commit(xy_xrows_type,ierr)
  call MPI_Type_commit(xy_yrows_type,ierr)

  call MPI_Type_commit(xz_xrows_type,ierr)
  call MPI_Type_commit(xz_zrows_type,ierr)

  call MPI_Type_commit(yz_yrows_type,ierr)
  call MPI_Type_commit(yz_zrows_type,ierr)
!-----------------------------------------------------------------------------------------
! filter ghost planes

  call MPI_Type_vector(ny*nz,    1+iforder/2 ,    nx, MPI_REAL8, fxrows_type, ierr)
  call MPI_Type_vector(   nz,    nx*(1+iforder/2), nx*ny, MPI_REAL8, fyrows_type, ierr)
  call MPI_Type_vector(    1, nx*ny*(1+iforder/2),     1, MPI_REAL8, fzrows_type, ierr)

  call MPI_Type_commit(fxrows_type,ierr)
  call MPI_Type_commit(fyrows_type,ierr)
  call MPI_Type_commit(fzrows_type,ierr)
!-----------------------------------------------------------------------------------------
! cshift end planes

  call MPI_Type_vector(ny*nz,     1,    nx, MPI_REAL8, xcs_type, ierr)
  call MPI_Type_vector(   nz,    nx, nx*ny, MPI_REAL8, ycs_type, ierr)
  call MPI_Type_vector(    1, nx*ny,     1, MPI_REAL8, zcs_type, ierr)

  call MPI_Type_commit(xcs_type,ierr)
  call MPI_Type_commit(ycs_type,ierr)
  call MPI_Type_commit(zcs_type,ierr)
!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    call write_header(io,'-')
  endif
!----------------------------------------------------------------------------------------
  return
  end subroutine initialize_topology
!=========================================================================================
  subroutine check_term_status(io,flag)
!=========================================================================================
! routine checks term_status and terminates cleanly if there is a problem
!-----------------------------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io    !io unit
  integer flag  !0 for simple termination, 1 for termination with savefiles
!-----------------------------------------------------------------------------------------
! broadcast status

  call MPI_Bcast(term_status,1,MPI_INTEGER,0,gcomm,ierr)

! terminate run for all processors

  if(term_status.eq.1) call terminate_run(io,flag)  !must be called by all processors
!----------------------------------------------------------------------------------------
  return
  end subroutine check_term_status
!-----------------------------------------------------------------------------------------
  end module topology_m

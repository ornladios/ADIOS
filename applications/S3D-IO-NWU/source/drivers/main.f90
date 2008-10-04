!========================================================================================
  program s3d
!========================================================================================
! Direct Numerical Simulation of the 3-D compressible Navier-Stokes
! equations with detailed chemical reactions
!----------------------------------------------------------------------------------------
  use param_m
  use topology_m
  use chemkin_m
  use runtime_m

  implicit none
!----------------------------------------------------------------------------------------
! integer variable declarations

  integer io          !output io unit
!----------------------------------------------------------------------------------------
! start initialization of MPI parameters

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
  call MPI_Comm_dup (MPI_COMM_WORLD, gcomm,ierr)
!----------------------------------------------------------------------------------------
! read data from s3d.in

  call read_input(6)
!----------------------------------------------------------------------------------------
! set output to screen or file

  if(i_write.eq.0) then
    io=6
  else
    io=7
    if(myid.eq.0) then
      open(unit=io,file='../data/'//trim(run_title)//'.out',status='REPLACE')
    endif
  endif
!----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    call write_header(io,'=')
    write(io,*) 'Hello, and welcome to S3D!'
    call write_header(io,'=')
  endif
!----------------------------------------------------------------------------------------
! get start date and time and write to file

  if(myid.eq.0) then
    call write_header(io,'-')
    call date_and_time(dat_1,tim_1)
    write(io,*) 'this run started at:'
    write(io,*)
    call write_date_and_time(dat_1,tim_1,io)
    call write_header(io,'-')
  endif
!----------------------------------------------------------------------------------------
! initialize primary modules
!----------------------------------------------------------------------------------------
! initialize chemkin module

  call initialize_chemkin(io,myid,ierr,gcomm)

! initialize param module (must be done after initializing chemkin module)

  call initialize_param(io,myid,ierr,gcomm)

! intialize topology module (must be done after initializing param module)

  call initialize_topology(io,nx,ny,nz,npx,npy,npz,iorder,iforder)

! sync processors

  call MPI_Barrier(gcomm,ierr)
!----------------------------------------------------------------------------------------
! start tasks here
!----------------------------------------------------------------------------------------
! mode of operation

  if(mode.eq.0) then          !solve governing equations
    call solve_driver(io)
  else
    call post_driver(io)      !post process results
  endif
!----------------------------------------------------------------------------------------
! get final date and time and write to file

  if(myid.eq.0) then
    call write_header(io,'-')
    write(io,*) 'S3D started at:'
    write(io,*)
    call write_date_and_time(dat_1,tim_1,io)
    write(io,*)
    write(io,*) 'S3D ended at:'
    write(io,*)
    call date_and_time(dat_2,tim_2)
    call write_date_and_time(dat_2,tim_2,io)
  endif
!----------------------------------------------------------------------------------------
! write closing statement

  if(myid.eq.0) then
    call write_header(io,'=')
    write(io,*) 'Program finished normally, goodbye!'
    call write_header(io,'=')
  endif

! close output file if it exists

  if(i_write.ne.0) then
    if(myid.eq.0) close(io)
  endif
!----------------------------------------------------------------------------------------
! deallocate arrays

  call allocate_chemkin_arrays(-1)
!----------------------------------------------------------------------------------------
! cleanup MPI processes

  call MPI_Finalize(ierr)
!----------------------------------------------------------------------------------------
  end program s3d

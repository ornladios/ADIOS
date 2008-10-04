!========================================================================================
  subroutine post_driver(io)
!========================================================================================
! routine for post-processing DNS data
!----------------------------------------------------------------------------------------
  use param_m                   ! initialized in main.f90
  use topology_m                ! initialized in main.f90
  use chemkin_m                 ! initialized in main.f90

  use runtime_m, only : i_time, time
  use variables_m
  use work_m
  use grid_m
  use reference_m
!  use transport_m
  use derivative_m
  use thermchem_m
  use runtime_m, only : run_title

  implicit none

!---------------------------- START variable declarations -------------------------------
  integer, intent(in) :: io         ! output io unit
  integer :: io_logfile = 15
  character*100 :: filename
  logical :: exists, file_end=.false.
!----------------------------- END variable declarations --------------------------------


!------------------------------- START initialization -----------------------------------
  call initialize_derivative(io)                ! initialize derivative module
  ! Ramanan - 07/16/06 Reference must be initialized before transport
  call initialize_reference(io,0.708)              ! must be done after transport
!  call initialize_transport(io)                 ! initialize transport module
  call initialize_grid(io)                      ! must be done after reference
  call initialize_thermchem(io)                 ! must be done after reference

  n_reg = 1

  call allocate_variables_arrays( 1 )           ! allocate memory for variables
  call allocate_work_arrays( 1 )                ! allocate memory for work arrays

!-------------------------------- END initialization ------------------------------------

! open the file containing the list of files to be post-processed
  if(myid==0) then
     filename = '../data/'//trim(run_title)//'.savefile.log'
     inquire(file=filename, exist=exists)
     if(.not. exists) then
        write(io,*)
        write(io,*)'ERROR: file does not exist:',trim(filename)
        write(io,*)'ABORTING!!!'
        term_status=1
     endif
  endif
  call check_term_status(io,0)

  if(myid==0) then
     open(unit=io_logfile, file=filename, status='old', form='formatted')
   ! chew up the header lines of this file
     read(io_logfile,*) 
     read(io_logfile,*) 
  endif
!--------------------------- START loop over all time steps -----------------------------
  do   ! keep doing this until the log file runs out of entries...

   ! myid=0 reads the log file to determine the next file to post-process

     if(myid==0) then

        read(io_logfile,50,end=60) filename, i_time, time
50      format(a50,4x,i7,8x,1pe9.3)

      ! if we got here, we successfully read a line from the log file
      ! so set the "file_end" flag to false

        file_end = .false.

        goto 30

      ! error trapping on read statement

60      file_end = .true.

30      continue

        if( file_end .or. (i_time<0) .or. (time<0) ) then

           write(io,*)
           call write_header( io, '*' )
           write(io,*) 'END of log file found.  Terminating post-processing run now.'
           call write_header( io, '*' )

        endif

      ! nondimensionalize the time
        time = time/time_ref

     endif
!---------------------------------------------------------------------
!    broadcast the file status - if the end of the file was reached, 
!    we will gracefully exit the code.

     call MPI_Bcast(file_end,1,MPI_LOGICAL,0,gcomm,ierr)
     if(file_end) then
        call terminate_run(io,0)
     endif
!---------------------------------------------------------------------
!    broadcast the information on the file to be post-processed 

     call MPI_Bcast(filename, 100, MPI_CHARACTER, 0, gcomm, ierr )
     call MPI_Bcast(i_time,     1, MPI_INTEGER,   0, gcomm, ierr )
     call MPI_Bcast(time,       1, MPI_REAL8,     0, gcomm, ierr )

!---------------------------------------------------------------------
!    read the field file for this timestep

     call read_savefile( io )  ! read the restart file containing the primitive variables
!     call set_variables( io )  ! set other variables...



!~~~~~~~~~~~~~~~~~~~ START user post processing routines ~~~~~~~~~~~~~~~~~~~

     !call write_tecplot_skip(io,1)
!     call write_basic_tecplot_file(io,0)

!~~~~~~~~~~~~~~~~~~~~ END user post processing routines ~~~~~~~~~~~~~~~~~~~~

  enddo
  if(myid==0) close( io_logfile )
!--------------------------- END of loop over all time steps ----------------------------


!------------------------------ START deallocate arrays ---------------------------------
  call allocate_variables_arrays(-1)
  call allocate_work_arrays(-1)
  call allocate_grid_arrays(-1)
!------------------------------- END deallocate arrays ----------------------------------
  return
  end subroutine post_driver

!=========================================================================================
  subroutine terminate_run(io,flag)
!=========================================================================================
! routine terminates run cleanly with respect to MPI and writing savefiles
! this routine MUST be called from outside any if(myid.eq.0) statements
! in other words, it should be called from all processors
!-----------------------------------------------------------------------------------------
  use topology_m, only : myid, ierr, gcomm
  use param_m, only : dat_1, dat_2, tim_1, tim_2

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io    !io unit
  integer flag  !0 for simple termination, 1 for termination with write savefile
!-----------------------------------------------------------------------------------------
! simple termination

  if(flag.eq.0) then

    if(myid.eq.0) write(io,*)

    call MPI_Barrier(gcomm,ierr)

    if(myid==0) then
      write(io,*) 'terminating run with a simple termination...'
      write(io,*)
    endif

    call MPI_Barrier(gcomm,ierr)

    if(myid.eq.0) then
      write(io,*) 'S3D started at:'
      write(io,*)
      call write_date_and_time(dat_1,tim_1,io)
      write(io,*)
      write(io,*) 'S3D ended at:'
      write(io,*)
      call date_and_time(dat_2,tim_2)
      call write_date_and_time(dat_2,tim_2,io)
    endif

    call MPI_Barrier(gcomm,ierr)

    if(myid.eq.0) call write_header(io,'=')

    call MPI_Finalize(ierr)

    stop

  endif
!-----------------------------------------------------------------------------------------
! termination with restart dump

  if(flag.eq.1) then

    if(myid.eq.0) write(io,*)

    call MPI_Barrier(gcomm,ierr)

    if(myid==0) then
      write(io,*) 'terminating run with a savefile...'
      write(io,*)
    endif

    call MPI_Barrier(gcomm,ierr)

    call write_savefile(io)

    call MPI_Barrier(gcomm,ierr)

    if(myid.eq.0) then
      write(io,*) 'S3D started at:'
      write(io,*)
      call write_date_and_time(dat_1,tim_1,io)
      write(io,*)
      write(io,*) 'S3D ended at:'
      write(io,*)
      call date_and_time(dat_2,tim_2)
      call write_date_and_time(dat_2,tim_2,io)
    endif

    call MPI_Barrier(gcomm,ierr)

    if(myid.eq.0) call write_header(io,'=')

    call MPI_Finalize(ierr)

    stop

  endif
!-----------------------------------------------------------------------------------------
  return
  end

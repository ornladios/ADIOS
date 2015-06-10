!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

program read_adios_f
  implicit none
  include "mpif.h"
  integer :: gcnt, vcnt, acnt, tstart, tstop, vrank, vtype, vtimed, ierr
  integer :: comm,i,j,k,l,m
  integer*8 :: fh, gh
  integer, dimension(2) :: dims, start, readsize
  character (len=100), dimension(10000) :: vnamelist
  character (len=100), dimension(10000) :: gnamelist
  real, allocatable,dimension(:,:)    :: particles
  integer :: mype, ierror,pe_size, size_per_proc
  integer :: numargs
  character(len=256) :: fname
  integer           :: time,iter
  real*8 :: start_time, end_time, total_time
  real*8 :: sz,gbs


#ifndef __GFORTRAN__
!#ifndef __GNUC__
    interface
         integer function iargc()
         end function iargc
    end interface
!#endif
#endif

  call MPI_Init (ierr)
  comm = MPI_COMM_WORLD
  call mpi_comm_rank(mpi_comm_world,mype,ierror)
  call mpi_comm_size(mpi_comm_world,pe_size,ierror)
  
  if (mype==0) then
     numargs = iargc()
     
     if (numargs < 1) then
        call MPI_Finalize (ierr)
        call exit(1)
     end if
     call getarg(1,fname)
  end if

  call MPI_BCAST(fname,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  call adios_fopen (fh, trim(fname)//char(0), comm, ierr)
  call adios_inq_file (fh,gcnt,vcnt,acnt,tstart,tstop,gnamelist,ierr) 
  call adios_gopen (fh, gh, gnamelist(1), ierr) 
  call adios_inq_group(gh, vcnt, vnamelist, ierr)
  if (mype==0) then
     do i=1,vcnt
        !qwrite(6,*) i,vnamelist(i)
     end do
  end if

  start(1:2)=0
  call adios_inq_var (gh, "ptrackedi"//char(0),vtype, vrank, vtimed, dims, ierr)
  size_per_proc = dims(2)/pe_size
  start(1)=0
  start(2)=mype*size_per_proc
  readsize(1)= dims(1)
  readsize(2)= size_per_proc 
  if (mype==0) then
      open(600,file='times.dat',ACCESS="APPEND")
  end if
  allocate(particles(readsize(1),readsize(2)))
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  start_time = MPI_WTIME()
  call adios_get_var (gh, "ptrackedi"//char(0), particles, start, readsize, tstart, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end_time = MPI_WTIME()
  total_time = end_time - start_time
  sz = 4.0 * dims(1)*dims(2)/1024.0/1024.0/1024.0
  gbs = sz/total_time
  if (mype==0) write(6,100) trim(fname),pe_size,sz,total_time,gbs
100 format(a18,1x,i4,2x,f12.3,f10.2,f12.3)

      if (mype==0) write(600,*) 'time: ',pe_size,end_time-start_time
     if (mype==0) close(600)
     if (mype==0) flush(600)
  deallocate(particles)
  call adios_gclose(gh, ierr)
  call adios_fclose(fh, ierr)
  call MPI_Finalize (ierr)
end program read_adios_f

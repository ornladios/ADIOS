!  modifications
!  early 2005: Ramanan Sankaran added capability for io from directories
!                               instead of tar-balls
!  17-02-2005: Evatt Hawkes.  Bug fix of new io capability - bug on restarts
!
!
!=========================================================================================
  subroutine read_input(io)
!=========================================================================================
! routine reads input file s3d.in for basic run parameters and broadcasts them
!-----------------------------------------------------------------------------------------
  use topology_m

  use param_m, only : mode
  use param_m, only : nx_g, ny_g, nz_g, npx, npy, npz
  use param_m, only : vary_in_x, vary_in_y, vary_in_z
  use param_m, only : periodic_x, periodic_y, periodic_z
  use param_m, only : iorder, iforder

  use runtime_m, only : i_restart, i_time_end, i_time_save, time_save_inc, time_save
  use runtime_m, only : run_title, i_write, i_time_mon, i_time_res, i_time_tec

  !use turbulence_m, only : i_turbulence

  use bc_m, only : nrf_x0, nrf_xl, nrf_y0, nrf_yl, nrf_z0, nrf_zl, relax_ct

  use grid_m, only : xmin, ymin, zmin, xmax, ymax, zmax

  !use transport_m, only : pr

  !use filter_m, only : i_time_fil

  use reference_m, only : g_ref, a_ref, t_o, rho_ref, lambda_ref, mach_no, re_real

  use thermchem_m, only : i_react

! wkliao:
  use mpi_io_m, only : io_method

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io

! local declarations

  character*100 filename
  logical exist
  integer idummy
  real rdummy
!-----------------------------------------------------------------------------------------
! set file name and inquire

  filename='../input/s3d.in'
  call inquire_about_input_file(filename,io)
!-----------------------------------------------------------------------------------------
! open and read file

  if(myid.eq.0) then
!-----------------------------------------------------------------------------------------
! open file

  open(unit=1,file=trim(filename),status='old')
!-----------------------------------------------------------------------------------------
! read mode

  read(1,*)
  read(1,*)
  read(1,*)

  read(1,*) mode
!-----------------------------------------------------------------------------------------
! read grid dimension parameters

  read(1,*)
  read(1,*)
  read(1,*)

  read(1,*) nx_g
  read(1,*) ny_g
  read(1,*) nz_g
  read(1,*) npx
  read(1,*) npy
  read(1,*) npz
!-----------------------------------------------------------------------------------------
! read run-time parameters

  read(1,*)
  read(1,*)
  read(1,*)

  read(1,*) i_write
  read(1,*) i_restart
  read(1,*) i_time_end
  read(1,*) i_time_save
  read(1,*) time_save_inc
!-----------------------------------------------------------------------------------------
! read geometry parameters

  read(1,*)
  read(1,*)
  read(1,*)

  read(1,*) run_title
  read(1,*) vary_in_x
  read(1,*) vary_in_y
  read(1,*) vary_in_z
  read(1,*) periodic_x
  read(1,*) periodic_y
  read(1,*) periodic_z
  read(1,*) idummy
  !read(1,*) i_turbulence
  read(1,*) nrf_x0
  read(1,*) nrf_xl
  read(1,*) nrf_y0
  read(1,*) nrf_yl
  read(1,*) nrf_z0
  read(1,*) nrf_zl
  read(1,*) relax_ct
!-----------------------------------------------------------------------------------------
! read physical parameters

  read(1,*)
  read(1,*)
  read(1,*)

  read(1,*) xmin
  read(1,*) ymin
  read(1,*) zmin
  read(1,*) xmax
  read(1,*) ymax
  read(1,*) zmax
  read(1,*) mach_no
  read(1,*) re_real
  read(1,*) rdummy
!-----------------------------------------------------------------------------------------
! read numerics switches parameters

  read(1,*)
  read(1,*)
  read(1,*)

  read(1,*) i_react
  read(1,*) iorder
  read(1,*) i_time_mon
  read(1,*) i_time_res
  read(1,*) i_time_tec
  read(1,*) iforder
  !read(1,*) i_time_fil
  read(1,*) idummy
!-----------------------------------------------------------------------------------------
! read required reference values

  read(1,*)
  read(1,*)
  read(1,*)

  read(1,*) g_ref
  read(1,*) a_ref
  read(1,*) t_o
  read(1,*) rho_ref
  read(1,*) lambda_ref
!-----------------------------------------------------------------------------------------
! wkliao: read parameters for MPI I/O

  read(1,*)                 ! 3 comment lines
  read(1,*)
  read(1,*)

  read(1,*) io_method

!-----------------------------------------------------------------------------------------
! end of read

  close(1)

  endif
!-----------------------------------------------------------------------------------------
! broadcast mode

  call MPI_Bcast(mode    , 1, MPI_INTEGER, 0, gcomm, ierr)
!-----------------------------------------------------------------------------------------
! broadcast grid parameters

  call MPI_Bcast(nx_g    , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(ny_g    , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(nz_g    , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(npx     , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(npy     , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(npz     , 1, MPI_INTEGER, 0, gcomm, ierr)
!-----------------------------------------------------------------------------------------
! broadcast run-time parameters

  call MPI_Bcast(i_write        , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(i_restart      , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(i_time_end     , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(i_time_save    , 1, MPI_INTEGER, 0, gcomm, ierr)
  call MPI_Bcast(time_save_inc  , 1, MPI_REAL8  , 0, gcomm, ierr)
!-----------------------------------------------------------------------------------------
! broadcast geometry parameters

  call MPI_Bcast(run_title   , 20, MPI_CHARACTER, 0, gcomm, ierr)
  call MPI_Bcast(vary_in_x   ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(vary_in_y   ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(vary_in_z   ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(periodic_x  ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(periodic_y  ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(periodic_z  ,  1, MPI_INTEGER  , 0, gcomm, ierr)
!  call MPI_Bcast(i_turbulence,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(nrf_x0      ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(nrf_xl      ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(nrf_y0      ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(nrf_yl      ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(nrf_z0      ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(nrf_zl      ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(relax_ct    ,  1, MPI_REAL8    , 0, gcomm, ierr)
!-----------------------------------------------------------------------------------------
! broadcast physical parameters

  call MPI_Bcast(xmin   , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(ymin   , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(zmin   , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(xmax   , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(ymax   , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(zmax   , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(mach_no, 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(re_real, 1, MPI_REAL8 , 0, gcomm, ierr)
  !call MPI_Bcast(pr     , 1, MPI_REAL8 , 0, gcomm, ierr)
!-----------------------------------------------------------------------------------------
! broadcast numerics parameters

  call MPI_Bcast(i_react      ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(iorder       ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(i_time_mon   ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(i_time_res   ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(i_time_tec   ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  call MPI_Bcast(iforder      ,  1, MPI_INTEGER  , 0, gcomm, ierr)
  !call MPI_Bcast(i_time_fil   ,  1, MPI_INTEGER  , 0, gcomm, ierr)
!-----------------------------------------------------------------------------------------
! broadcast required reference values

  call MPI_Bcast(g_ref  , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(a_ref  , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(t_o    , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(rho_ref    , 1, MPI_REAL8 , 0, gcomm, ierr)
  call MPI_Bcast(lambda_ref , 1, MPI_REAL8 , 0, gcomm, ierr)
!-----------------------------------------------------------------------------------------
! wkliao: broadcast variables for I/ method to be usedO
  call MPI_Bcast(io_method,   1, MPI_INTEGER, 0, gcomm, ierr)

!-----------------------------------------------------------------------------------------
! set time_save equal to time_save_inc initially

  time_save = time_save_inc

  return
  end subroutine read_input
!=========================================================================================
  subroutine write_savefile(io)
!=========================================================================================
! writes data file for post-processing and restarting
!-----------------------------------------------------------------------------------------
  use topology_m, only : myid, gcomm, ierr
  use runtime_m, only : i_time, i_time_end, run_title, time, i_restart, time_restart, i_time_save
  use reference_m, only : time_ref

! wkliao: paramaters related to MPI-IO, pnetcdf, and hdf5
  use mpi_io_m,  only : io_method, mpi_file_io
  use pnetcdf_m, only : pnetcdf_write
  use hdf5_m,    only : hdf5_write
  use adios_m,   only : write_adios

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io

! other declarations

  character*5 myid_ext
  character*9 :: time_ext
  character*100  :: filename, tarcmd, tartgt, tarsrc, dirname
  logical :: exist
!-----------------------------------------------------------------------------------------

! BUG FIX Evatt Hawkes NOV 2004
! Barrier here should not be necessary, but recommended by PNNL consultants
! after we experienced some io issues on large runs
  call MPI_Barrier( gcomm, ierr )

  write(time_ext,'(1pe9.3)') time*time_ref       

! set file extention strings

!  call set_file_extension_string(myid_ext,myid,5)
! Change by Ramanan. The above subroutine does not do its job right
! on the Cray. Anyway... i dont think it needs a separate routine
! to format an integer into a string.. It should be as simple as this...
  write(myid_ext, '(I5.5)') myid

! wkliao: set a different file name for MPI-IO, pnetcdf, and HDF5 methods
! wkliao: no subdirectories are created if io_method > 0
if (io_method .EQ. 1) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.mpi'
elseif (io_method .EQ. 2) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.nc'
elseif (io_method .EQ. 3) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.h5'
elseif (io_method .EQ. 4) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.bp'
else

! set file name

#ifdef SAVEFILEINSEPDIR
  dirname = '../data/'//trim(run_title)//'.'//time_ext//'/'
  if(myid == 0) then
#ifdef SYSTEMCALLWONTWORK
! when system calls do not work, as on infiniband
      call makedirectory(trim(dirname)//char(0))
#else
! when system calls work ok
      tarcmd = 'mkdir -p '//trim(dirname)
      call execute_command( trim(tarcmd) )
#endif
  end if
  !All processes need to wait for the directory to be created
  call MPI_Barrier( gcomm, ierr )
  filename=trim(dirname)//'field.'//myid_ext
#else /*Not SAVEFILEINSEPDIR*/
  filename='field.'//myid_ext
#endif /*SAVEFILEINSEPDIR*/
endif

  if(myid.eq.0) then
    write(io,100) 'writing save    files for: i_time = ',i_time,  &
                  ', time = ',time*time_ref,' (sec)'
  endif

! wkliao:
if (io_method .GT. 0) then
    ! first delete the file if exists
    if (myid == 0) then
        inquire(file=trim(filename),exist=exist)
        if (exist) then
#ifdef SYSTEMCALLWONTWORK
            call unlink_file(trim(filename)//char(0))
#else
            call execute_command('rm -f '//trim(filename))
#endif
        endif
    endif
    call MPI_Barrier(gcomm, ierr)

    if (io_method .EQ. 1) then
        call mpi_file_io(filename, 'w')
    elseif (io_method .EQ. 2) then
        call pnetcdf_write(filename)
    elseif (io_method .EQ. 3) then
        call hdf5_write(filename)
    elseif (io_method .EQ. 4) then
        call write_adios(filename)
    endif
else
  open(unit=1,file=trim(filename),status='unknown',form='unformatted')

  call readwrite_savefile_data(io,1,'w')

  close(1)
endif
!-----------------------------------------------------------------------------------------
! file compression

  call MPI_Barrier( gcomm, ierr )

#ifndef SAVEFILEINSEPDIR

! set system commands for UNIX platforms

  tartgt = '../data/'//trim(run_title)//'.'//time_ext//'.tar'
  tarsrc = 'field.*'

  tarcmd = 'tar -cf '//trim(tartgt)//' '//trim(tarsrc)
  if(myid==0) call execute_command( trim(tarcmd) )
  endif
#ifdef SYSTEMCALLWONTWORK
  if(myid==0) call unlink_file(trim(tarsrc)//char(0))
#else
  if(myid==0) call execute_command('rm -f '//trim(tarsrc) )
#endif

! set system commands for PC platform (DOS)
! note: you must have a tar package installed and it must 
! reside in a directory contained in the DOS path
! recommended files are tar.exe and cygwin1.dll 

# if PC

    tartgt = '..\data\'//trim(run_title)//'.'//time_ext//'.tar' 
    tarsrc = 'field.*'

    if (myid == 0) then
      tarcmd = 'tar -cf '//trim(tartgt)//' '//trim(tarsrc)
      call execute_command( trim(tarcmd) )      ! tar the files
      tarcmd = 'del '//trim(tarsrc)
      call execute_command( trim(tarcmd) )      ! delete the old files
    endif

# endif

#endif /* SAVEFILEINSEPDIR */

! write header

  if(myid.eq.0) then
    call write_header(io,'-')
  endif
!-----------------------------------------------------------------------------------------
! write log file to keep track of when save files were written

  if(myid==0) then

    filename='../data/'//trim(run_title)//'.savefile.log'
    inquire(file=trim(filename),exist=exist)

    if(exist) then

      open(unit=1,file=trim(filename),status='old',position='append')

    else

      open(unit=1,file=trim(filename),status='unknown')
      write(1,51) 'file name','time step','time (sec)'
      write(1,52)

    endif

! wkliao:
if (io_method .EQ. 1) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.mpi'
elseif (io_method .EQ. 2) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.nc'
elseif (io_method .EQ. 3) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.h5'
elseif (io_method .EQ. 4) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.bp'
else
    filename = trim(run_title)//'.'//time_ext//'.tar'
endif

    write(1,50) filename, i_time, time*time_ref

    close(1)

  endif
!-----------------------------------------------------------------------------------------
! format statements

  50 format(a50,4x,i7,8x,1pe9.3)
  51 format(a25,27x,a9,7x,a10)
  52 format(70('-'))
  100 format(1x,a,i7,a,1pe9.3,a)
!-----------------------------------------------------------------------------------------
  return
  end subroutine write_savefile
!=========================================================================================
  subroutine read_savefile(io)
!=========================================================================================
! reads data file for post-processing and restarting
!-----------------------------------------------------------------------------------------
  use topology_m, only : myid, gcomm, ierr
  use runtime_m, only : i_time, run_title, time, i_restart
  use param_m, only : mode
  use reference_m, only : time_ref

! wkliao: paramaters related to MPI-IO, pnetcdf, and hdf5
  use mpi_io_m,  only : io_method, mpi_file_io
  use pnetcdf_m, only : pnetcdf_read
  use hdf5_m,    only : hdf5_read


  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io

! other declarations

  character*5 myid_ext
  character*9 :: time_ext
  character*100 :: filename, tartgt, tarcmd
  logical exist
!-----------------------------------------------------------------------------------------
! set file extention strings

! call set_file_extension_string(myid_ext,myid,5)
! Change by Ramanan. The above subroutine does not do its job right
! on the Cray. Anyway... i dont think it needs a separate routine
! to format an integer into a string.. It should be as simple as this...
  write(myid_ext, '(I5.5)') myid

  write(time_ext,'(1pe9.3)') time*time_ref       !set time stamp
!-----------------------------------------------------------------------------------------
! file decompression
! only untar files on post-processing runs (mode=1)
! on restarts, the user should untar the files before starting the code

#ifndef SAVEFILEINSEPDIR

  if(mode==1) then     


!   system commands for PC (DOS) platform

#   if PC

      tartgt = '..\data\'//trim(run_title)//'.'//time_ext//'.tar' 
      tarcmd = 'tar -xf '//trim(tartgt)
      if(myid==0)  call execute_command( trim(tarcmd) )

#   endif

!   system commands for UNIX platforms

    tartgt = '../data/'//trim(run_title)//'.'//time_ext//'.tar'
    tarcmd = 'tar -xf '//trim(tartgt)
    if(myid==0)  call execute_command( trim(tarcmd) )

  endif
#endif /*SAVEFILEINSEPDIR*/

  call MPI_Barrier( gcomm, ierr )
!-----------------------------------------------------------------------------------------
! set file name

! wkliao: MPI-IO, pnetcdf, and hdf5 methods
if (io_method .EQ. 1) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.mpi'
elseif (io_method .EQ. 2) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.nc'
elseif (io_method .EQ. 3) then
    filename = '../data/'//trim(run_title)//'.'//time_ext//'.field.h5'
else
  if(i_restart==1 .and. mode==0) then
     filename='../data/'//'field.'//myid_ext
  else
#ifdef SAVEFILEINSEPDIR
     filename='../data/'//trim(run_title)//'.'//time_ext//'/field.'//myid_ext
#else
     filename='field.'//myid_ext
#endif /*SAVEFILEINSEPDIR*/
  endif
endif

! inquire about file existence

  inquire(file=trim(filename),exist=exist)

  if(.not.exist) then
!jcs added capability to restart from last written file... (2-14-03)
     if(i_restart==1) then
        call restartFromLastFile(io)
#ifdef SAVEFILEINSEPDIR
        !bug fix Evatt Hawkes 17-FEB-2005
        write(time_ext,'(1pe9.3)') time*time_ref   
! wkliao:
!         filename='../data/'//trim(run_title)//'.'//time_ext//'/field.'//myid_ext
! #else
!         filename='field.'//myid_ext
#endif
     else
        call MPI_Barrier(gcomm,ierr)
        write(io,*) 'the following restart file does not exist'
        write(io,*) trim(filename)
        call MPI_Barrier(gcomm,ierr)
        call terminate_run(io,0)
     endif
  else
    if(myid.eq.0) then
      write(io,*) 'reading save files...'
      write(io,*)
    endif
  endif

  call MPI_Barrier(gcomm,ierr)

! wkliao: MPI I/O
if (io_method .GT. 0) then
    if (io_method .EQ. 1) then
        call mpi_file_io(filename, 'r')
    elseif (io_method .EQ. 2) then
        call pnetcdf_read(filename)
    elseif (io_method .EQ. 3) then
        call hdf5_read(filename)
    endif
else
! open file

!  open(unit=1,file=trim(filename),convert='BIG_ENDIAN',status='old',form='unformatted')
  open(unit=1,file=trim(filename),status='old',form='unformatted')

! read data

  call readwrite_savefile_data(io,1,'r')

! close file

  close(1)
endif
!-----------------------------------------------------------------------------------------
! delete untarred field files

  call MPI_Barrier(gcomm,ierr)

#ifndef SAVEFILEINSEPDIR

#if PC

    tarcmd = 'del '//trim(filename)

#else

     tarcmd = 'rm '//trim(filename)

#endif

#ifdef SYSTEMCALLWONTWORK
     call unlink_file(trim(filename)//char(0))
#else
     call execute_command( trim(tarcmd) )
#endif

#endif /*SAVEFILEINSEPDIR*/
!-----------------------------------------------------------------------------------------
! write header

  if(myid.eq.0) then
    call write_header(io,'-')
  endif
!-----------------------------------------------------------------------------------------
  return
  end subroutine read_savefile
!=========================================================================================
  subroutine restartFromLastFile(io)
!=========================================================================================
    ! restarts the code from the last file written
    ! jcs 2-14-03
    use topology_m
    use runtime_m, only : i_time, time, run_title
    use reference_m, only : time_ref
    implicit none
    integer, intent(in) :: io
    real :: time_old
    integer, parameter :: io_logfile = 20
    character*100 :: filename, tarcmd, tartgt
    character*9 :: time_ext
    logical :: exists

    filename = '../data/'//trim(run_title)//'.savefile.log'
    call inquire_about_input_file(filename,io)
    if(myid==0) then
       open(unit=io_logfile, file=filename, status='old', form='formatted')
       ! chew up the header lines of this file
       read(io_logfile,*)
       read(io_logfile,*)
       do
          read(io_logfile,50,end=60) filename, i_time, time
          write(*,*)trim(filename),i_time,time
          time_old = time
       enddo
50     format(a50,4x,i7,8x,1pe9.3)
60     continue
       close(io_logfile)
    endif
    call MPI_Bcast(filename, 100, MPI_CHARACTER, 0, gcomm, ierr )
    call MPI_Bcast(i_time,     1, MPI_INTEGER,   0, gcomm, ierr )
    call MPI_Bcast(time,       1, MPI_REAL8,     0, gcomm, ierr )

  !-- nondimensionalize the time
    time = time/time_ref
    write(time_ext,'(1pe9.3)') time*time_ref       !set time stamp for tar files

  !-- untar the last file found in savefile.log
!bug  fix Evatt Hawkes 17-FEB-2005
#ifndef SAVEFILEINSEPDIR

# if PC

    tartgt = '..\data\'//trim(run_title)//'.'//time_ext//'.tar' 
    tarcmd = 'tar -xf '//trim(tartgt)
    if(myid==0)  call execute_command( trim(tarcmd) )

# else

!   system commands for UNIX platforms
    tartgt = '../data/'//trim(run_title)//'.'//time_ext//'.tar'
    tarcmd = 'tar -xf '//trim(tartgt)
    if(myid==0)  call execute_command( trim(tarcmd) )
#endif

#else
#endif

    call MPI_Barrier(gcomm,ierr)
    return
  end subroutine restartFromLastFile
!=========================================================================================
  subroutine readwrite_savefile_data(io,io_savefile,input)
!=========================================================================================
! reads and writes specific data for save file
! this routine exists so that the reading and writing are both in one location
!-----------------------------------------------------------------------------------------
  use topology_m, only : myid

  use runtime_m, only : time, tstep, time_save
  use variables_m, only : temp, pressure, yspecies, u
  use bc_m, only : qx_bc, qy_bc, qz_bc, pout
  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io, io_savefile
  character*1 input
!-----------------------------------------------------------------------------------------
! write data

  if(input.eq.'w') then

    write(io_savefile) time
    write(io_savefile) tstep
    write(io_savefile) time_save
    write(io_savefile) yspecies
    write(io_savefile) temp
    write(io_savefile) pressure
    write(io_savefile) u
    write(io_savefile) pout

  elseif(input.eq.'r') then

    read(io_savefile) time
    read(io_savefile) tstep
    read(io_savefile) time_save
    read(io_savefile) yspecies
    read(io_savefile) temp
    read(io_savefile) pressure
    read(io_savefile) u
    read(io_savefile) pout

  else

    if(myid.eq.0) then
      write(io,*) 'improper setting for variable input'
      write(io,*) 'in routine readwrite_savefile_data'
    endif
    call terminate_run(io,0)  !must be called by all processors

  endif
!-----------------------------------------------------------------------------------------
  return
  end subroutine readwrite_savefile_data
!=========================================================================================
  subroutine write_header(io,char_in)
!=========================================================================================
! writes header for neatness
!-----------------------------------------------------------------------------------------
  implicit none

! declarations passed in

  integer io
  character*1 char_in

! other declarations

  integer num
  parameter(num=80)
  character*1 char(num)
  integer i
!-----------------------------------------------------------------------------------------
! fill char

  char(:)=char_in

! write header

  write(io,1) ' ', (char(i), i=1,num,1)
1 format(a1,100(a1))
!-----------------------------------------------------------------------------------------
  return
  end subroutine write_header
!=========================================================================================
  subroutine inquire_about_input_file(filename,io)
!=========================================================================================
! routine inquires about the existence of an input file
! and terminates the run cleanly if it doesn't exist
! this routine MUST be called from outside of any if(myid.eq.0) statements
! routine works best when filename passed in character*100 (which is the convention)
!-----------------------------------------------------------------------------------------
  use topology_m

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  character*100 filename
  integer io

! local declarations

  logical exist
!-----------------------------------------------------------------------------------------
! inquire about existence of file for myid=0

  if(myid.eq.0) then

    inquire(file=trim(filename),exist=exist)

    if(.not.exist) then   !does not exist

      write(io,*) 'the following input file does not exist:'
      write(io,*) trim(filename)

      term_status=1

    endif

  endif

! check termination status and terminate if file doesn't exist

  call check_term_status(io,0)
!----------------------------------------------------------------------------------------
  return
  end subroutine inquire_about_input_file
!=========================================================================================
  subroutine set_file_extension_string(file_ext,num,max)
!=========================================================================================
! converts a positive integer (num) to character (file_ext) of length max
! in other words, 0 < num < (10**max)-1
!
! routine uses function getchar (see below)
!
! note: this routine takes advantage of the fortran convention for passing arrays
!   through argument lists; in particular, in this routine file_ext is delcared
!   as a character array of length equal to the length of the single character
!   declaration of the calling routine; i'm sorry to have to do this, but it
!   made things very nice in this routine; in short:
!
!   character*max file_ext => character file_ext(max)
!-----------------------------------------------------------------------------------------
  implicit none

! declarations passed in

  character(*) file_ext
  integer num, max

! other declarations

  integer i, m, n, temp
  character getchar
!-----------------------------------------------------------------------------------------
! check for postiveness

  if(num.lt.0) then
    write(6,*) 'num passed into routine'
    write(6,*) 'set_file_extension_string is negative'
    call terminate_run(6,0)
  endif

! check for maximum value of num

  if(num.gt.(10**max)-1) then
    write(6,*) 'num passed into routine set_file_extension_string'
    write(6,*) 'is greater than the character file_ext will allow'
    call terminate_run(6,0)
  endif
!-----------------------------------------------------------------------------------------
! zero stuff (must do this!)

  n=0; temp=0

  do i=1,max,1
    file_ext(i:i)='0'
  enddo

! set file extension

  do i=max,1,-1
    temp=n
    n=(num/10**(i-1))
    m=max-i+1
    file_ext(m:m)=getchar(int(n-(temp*10)))
  enddo
!-----------------------------------------------------------------------------------------
  return
  end subroutine set_file_extension_string
!=========================================================================================
  function getchar(n)
!=========================================================================================
! getchar returns a character corresponding to n
!-----------------------------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  character getchar
  integer n
!-----------------------------------------------------------------------------------------
  if(n.eq.0) getchar = '0'
  if(n.eq.1) getchar = '1'
  if(n.eq.2) getchar = '2'
  if(n.eq.3) getchar = '3'
  if(n.eq.4) getchar = '4'
  if(n.eq.5) getchar = '5'
  if(n.eq.6) getchar = '6'
  if(n.eq.7) getchar = '7'
  if(n.eq.8) getchar = '8'
  if(n.eq.9) getchar = '9'
!-----------------------------------------------------------------------------------------
  return
  end function getchar
!=========================================================================================
  subroutine write_date_and_time(dat,tim,io)
!=========================================================================================
! writes date and time in formatted sequence to file io
!-----------------------------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------------------------
! declarations

  integer io

  character*1 dat(8), tim(10)
!-----------------------------------------------------------------------------------------
! write time and date

  write(io,1) dat(5),dat(6),dat(7),dat(8),dat(1),dat(2),dat(3),dat(4)
  write(io,2) tim(1),tim(2),tim(3),tim(4),tim(5),tim(6)

  1 format(1x,'date = ',a1,a1,'/',a1,a1,'/',a1,a1,a1,a1)
  2 format(1x,'time = ',a1,a1,':',a1,a1,':',a1,a1)
!-----------------------------------------------------------------------------------------
  return
  end subroutine write_date_and_time
!=========================================================================================
  subroutine tar_input_files(io)
!=========================================================================================
! routine tars input files after initialization for record of run parameters
! places tar file in data directory
!
! PC users must have tar.exe and cygwin1.dll files in path directory
!-----------------------------------------------------------------------------------------
  use topology_m, only : myid
  use runtime_m, only : run_title

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io

! local declarations

  character*100 tarcmd, copycmd1, copycmd2, delcmd
!-----------------------------------------------------------------------------------------
  if(myid == 0) then

    write(io,*) 'tarring input files for record of run parameters...'
    write(io,*) 'tar file is located at ../data/'//trim(run_title)//'.in.tar'
    write(io,*) 

#   if PC
      copycmd1 = 'copy ..\input\chem.asc ..\run'
      copycmd2 = 'copy ..\input\*.in ..\run'
      tarcmd = 'tar -cf ..\data\'//trim(run_title)//'.in.tar *.in chem.asc'
      delcmd = 'del *.in chem.asc'
#   else
      copycmd1 = 'cp ../input/chem.asc ../run'
      copycmd2 = 'cp ../input/*.in ../run'
      tarcmd = 'tar -cf ../data/'//trim(run_title)//'.in.tar *.in chem.asc'
      delcmd = 'rm *.in chem.asc'
#   endif

      call execute_command( trim(copycmd1) )
      call execute_command( trim(copycmd2) )
      call execute_command( trim(tarcmd) )
      call execute_command( trim(delcmd) )

    call write_header(io,'-')

  endif
!-----------------------------------------------------------------------------------------
  return
  end subroutine tar_input_files
!=========================================================================================
  subroutine generate_active_file(io)
!=========================================================================================
! routine generates active input file for changing a few key parameters
! while the code is running
!-----------------------------------------------------------------------------------------
  use topology_m, only : myid, term_status
  use runtime_m, only : run_title, i_time_end, i_time_save, time_save_inc
  use runtime_m, only : i_time_mon, i_time_res, i_time_tec
  !use filter_m, only : i_time_fil
  use rk_m, only : cont_switch, cfl_switch, tstep_min, tstep_max, tstep_init
  use rk_m, only : rk_rtol, rk_atol
  use rk_m, only : i_time_cont

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io

! local declarations

  character*100 filename
!-----------------------------------------------------------------------------------------
! generate active.in file

  if(myid.eq.0) then

!   write header

    write(io,*) 'generating active.in file...'
    write(io,*)

!   open file

    filename='../input/active.in'
    open(unit=56,file=trim(filename),status='unknown')

!   begin general parameters

    write(56,1) '0-continue, 1-terminate soon with savefiles         (term_status) ',  &
                term_status
    write(56,1) 'ending time step                                     (i_time_end) ',  &
                i_time_end
    write(56,1) 'frequency to save fields in restart files           (i_time_save) ',  &
                i_time_save
    write(56,2) 'time period to save fields in restart files       (time_save_inc) ',  &
                time_save_inc
    write(56,1) 'frequency to monitor min/max files and active        (i_time_mon) ',  &
                i_time_mon
    write(56,1) 'frequency to check resolution; set < 0 for no check  (i_time_res) ',  &
                i_time_res
    write(56,1) 'frequency write current tecplot file                 (i_time_tec) ',  &
                i_time_tec
    !write(56,1) 'frequency to filter solution vector                  (i_time_fil) ',  &
    !            i_time_fil

!   begin controller parameters

    if(cont_switch==0) then

    write(56,2) 'initial/constant timestep (sec)                      (tstep_init) ',  &
                tstep_init

    else

    write(56,2) 'minimum timestep for controller (sec)                 (tstep_min) ',  &
                tstep_min
    write(56,2) 'maximum timestep for controller (sec)                 (tstep_max) ',  &
                tstep_max
    write(56,2) 'relative Runge-Kutta error tolerance                    (rk_rtol) ',  &
                rk_rtol
    write(56,2) 'absolute Runge-Kutta error tolerance                    (rk_atol) ',  &
                rk_atol
    write(56,1) 'controller cfl check, 0=off, 1=on                    (cfl_switch) ',  &
                cfl_switch

    endif

    write(56,1) 'timestep frequency to write controller info         (i_time_cont) ',  &
                i_time_cont

!   close file

    close(56)

!   write header

    call write_header(io,'-')

  endif
!----------------------------------------------------------------------------------------
! format statements

  1 format(a66,i12)
  2 format(a66,1pe12.5)
!----------------------------------------------------------------------------------------
  return
  end subroutine generate_active_file
!=========================================================================================
  subroutine read_active_file(io)
!=========================================================================================
! routine reads active input file and changes key parameters
! while the code is running
!
!
! BUG FIX Evatt Hawkes
! The previous code (whole subroutine commented out below) did not compile on seaborg.
! The new code is less efficient but more portable.
!
!-----------------------------------------------------------------------------------------
  use topology_m
  use runtime_m, only : run_title, i_time_end, i_time_save, time_save_inc
  use runtime_m, only : i_time_mon, i_time_res, i_time, i_time_tec
  !use filter_m, only : i_time_fil
  use rk_m, only : cont_switch, cfl_switch, tstep_min, tstep_max, tstep_init
  use rk_m, only : rk_rtol, rk_atol
  use rk_m, only : i_time_cont

  implicit none
!-----------------------------------------------------------------------------------------
! declarations passed in

  integer io

! local declarations

  character*100 filename
  logical exist
  character*66 dummy  !length specific to format statement in routine generate_active_file

  integer, parameter :: n_intgrs=9, n_reals=6

  real, dimension(n_reals) :: r
  integer, dimension(n_intgrs) :: i

!-----------------------------------------------------------------------------------------
! return if not desired

  if(mod(i_time,i_time_mon).ne.0) return
!-----------------------------------------------------------------------------------------
! check for existence of active.in file

  if(myid.eq.0) then

    filename='../input/active.in'
    inquire(file=trim(filename),exist=exist)

    if(exist) then   !continue with read

      continue

    else      !generate and return

      write(io,*) 'active.in file does not exist...'
      write(io,*) 'regenerating active.in file...'
      call generate_active_file(io)
      return

    endif

  endif
!-----------------------------------------------------------------------------------------
! read active.in file

!jcs
!!$  time_save_inc=0.0
!!$  i_time_save=0

  if(myid.eq.0) then

!   open file

    filename='../input/active.in'
    open(unit=56,file=trim(filename),status='unknown',err=10)

!   begin general parameters

    read(56,1) dummy, term_status
    read(56,1) dummy, i_time_end
    read(56,1) dummy, i_time_save
    read(56,2) dummy, time_save_inc
    read(56,1) dummy, i_time_mon
    read(56,1) dummy, i_time_res
    read(56,1) dummy, i_time_tec
    !read(56,1) dummy, i_time_fil
    read(56,1) dummy, i_time_tec

!   begin controller parameters

    if(cont_switch==0) then
      read(56,2) dummy, tstep_init
    else
      read(56,2) dummy, tstep_min
      read(56,2) dummy, tstep_max
      read(56,2) dummy, rk_rtol
      read(56,2) dummy, rk_atol
      read(56,1) dummy, cfl_switch
    endif

    read(56,1) dummy, i_time_cont
    goto 20

10  write(io,*) 'there is an error reading the active.in file...'
    write(io,*) 'proceeding without update of active variables...'

20  continue

    close(56)

  endif
!----------------------------------------------------------------------------------------
! broadcast parameters

  i(1) = term_status;    r(1) = time_save_inc
  i(2) = i_time_end ;    r(2) = tstep_init
  i(3) = i_time_save;    r(3) = tstep_min
  i(4) = i_time_mon ;    r(4) = tstep_max
  i(5) = i_time_res ;    r(5) = rk_rtol
  i(6) = i_time_tec ;    r(6) = rk_atol
  !i(7) = i_time_fil
  i(8) = cfl_switch
  i(9) = i_time_cont

  call MPI_Bcast( i, n_intgrs, MPI_INTEGER, 0, gcomm, ierr )
  call MPI_Bcast( r, n_reals,  MPI_REAL8, 0, gcomm, ierr )

  term_status = i(1);     time_save_inc = r(1)
  i_time_end  = i(2);     tstep_init    = r(2)
  i_time_save = i(3);     tstep_min     = r(3)
  i_time_mon  = i(4);     tstep_max     = r(4)
  i_time_res  = i(5);     rk_rtol       = r(5)
  i_time_tec  = i(6);     rk_atol       = r(6)
  !i_time_fil  = i(7);
  cfl_switch  = i(8);
  i_time_cont = i(9);
!----------------------------------------------------------------------------------------
! terminate if desired

  if(term_status.eq.1) call terminate_run(io,1)
  return
!----------------------------------------------------------------------------------------
! format statements

  1 format(a66,i12)
  2 format(a66,1pe12.5)
!----------------------------------------------------------------------------------------
  return
  end subroutine read_active_file


! BUG FIX Evatt Hawkes
! The following code is theoretically more efficient since only 
! one broadcast is made, but it does not compile on some platforms (seaborg).
!
!!=========================================================================================
!  subroutine read_active_file(io)
!!=========================================================================================
!! routine reads active input file and changes key parameters
!! while the code is running
!!
!! note: there sometimes is an address error on the last real being broadcast
!! in the datatype and the last real doesn't get broadcast properly
!! thus the last real is always a dummy real to avoid this problem
!!-----------------------------------------------------------------------------------------
!  use topology_m
!  use runtime_m, only : run_title, i_time_end, i_time_save, time_save_inc
!  use runtime_m, only : i_time_mon, i_time_res, i_time, i_time_tec
!  use filter_m, only : i_time_fil
!  use rk_m, only : cont_switch, cfl_switch, tstep_min, tstep_max, tstep_init
!  use rk_m, only : rk_rtol, rk_atol
!  use rk_m, only : i_time_cont
!
!  implicit none
!!-----------------------------------------------------------------------------------------
!! declarations passed in
!
!  integer io
!
!! local declarations
!
!  character*100 filename
!  logical exist
!  character*66 dummy  !length specific to format statement in routine generate_active_file
!
!  integer, parameter :: n_intgrs=9, n_reals=6
!
!  type passDataType
!     real, dimension(n_reals) :: r
!     integer, dimension(n_intgrs) :: i
!  end type passDataType
!
!  type(passDataType) :: psData
!  integer, save :: active_datatype
!  logical, save :: initialized_MPI_datatype = .false.
!
!  integer :: type_count(2), types(2)
!  integer :: displace(2)
!  integer :: displace_i
!!-----------------------------------------------------------------------------------------
!! return if not desired
!
!  if(mod(i_time,i_time_mon).ne.0) return
!!-----------------------------------------------------------------------------------------
!! check for existence of active.in file
!
!  if(myid.eq.0) then
!
!    filename='../input/active.in'
!    inquire(file=trim(filename),exist=exist)
!
!    if(exist) then   !continue with read
!
!      continue
!
!    else      !generate and return
!
!      write(io,*) 'active.in file does not exist...'
!      write(io,*) 'regenerating active.in file...'
!      call generate_active_file(io)
!      return
!
!    endif
!
!  endif
!!-----------------------------------------------------------------------------------------
!! read active.in file
!
!!jcs
!!!$  time_save_inc=0.0
!!!$  i_time_save=0
!
!  if(myid.eq.0) then
!
!   open file
!
!    filename='../input/active.in'
!    open(unit=56,file=trim(filename),status='unknown',err=10)
!
!!   begin general parameters
!
!    read(56,1) dummy, term_status
!    read(56,1) dummy, i_time_end
!    read(56,1) dummy, i_time_save
!    read(56,2) dummy, time_save_inc
!    read(56,1) dummy, i_time_mon
!    read(56,1) dummy, i_time_res
!    read(56,1) dummy, i_time_tec
!    read(56,1) dummy, i_time_fil
!
!!   begin controller parameters
!
!    if(cont_switch==0) then
!      read(56,2) dummy, tstep_init
!    else
!      read(56,2) dummy, tstep_min
!      read(56,2) dummy, tstep_max
!      read(56,2) dummy, rk_rtol
!      read(56,2) dummy, rk_atol
!      read(56,1) dummy, cfl_switch
!    endif
!
!    read(56,1) dummy, i_time_cont
!    goto 20
!
!10  write(io,*) 'there is an error reading the active.in file...'
!    write(io,*) 'proceeding without update of active variables...'
!
!20  continue
!
!    close(56)
!
!  endif
!!----------------------------------------------------------------------------------------
!! broadcast parameters
!
!  psData%i(1) = term_status;    psData%r(1) = time_save_inc
!  psData%i(2) = i_time_end ;    psData%r(2) = tstep_init
!  psData%i(3) = i_time_save;    psData%r(3) = tstep_min
!  psData%i(4) = i_time_mon ;    psData%r(4) = tstep_max
!  psData%i(5) = i_time_res ;    psData%r(5) = rk_rtol
!  psData%i(6) = i_time_tec ;    psData%r(6) = rk_atol
!  psData%i(7) = i_time_fil
!  psData%i(8) = cfl_switch
!  psData%i(9) = i_time_cont
!
!  ! create MPI datatype to allow us to pass all of these variables at once...
!  ! this should save some time
!  if (.not. initialized_MPI_datatype) then
!     type_count = (/ n_intgrs, n_reals /) !number of each type (in order)
!     types = (/ MPI_INTEGER, MPI_REAL8 /) !list of datatypes
!
!   ! compute memory displacements for datatypes in structure
!     call MPI_Type_extent( MPI_INTEGER, displace_i, ierr )
!
!   ! set displacements for each datatype in the structure
!     displace(1) = 0           
!     displace(2) = displace(1) + type_count(1)*displace_i
!
!   ! create and commit the MPI data structures
!     call MPI_Type_struct( 2, type_count, displace, types, active_datatype, ierr )
!     call MPI_Type_commit( active_datatype, ierr )
!
!     initialized_MPI_datatype = .true.
!  endif
!
!  call MPI_Bcast( psData, 1, active_datatype, 0, gcomm, ierr )
!
!  term_status = psData%i(1);     time_save_inc = psData%r(1)
!  i_time_end  = psData%i(2);     tstep_init    = psData%r(2)
!  i_time_save = psData%i(3);     tstep_min     = psData%r(3)
!  i_time_mon  = psData%i(4);     tstep_max     = psData%r(4)
!  i_time_res  = psData%i(5);     rk_rtol       = psData%r(5)
!  i_time_tec  = psData%i(6);     rk_atol       = psData%r(6)
!  i_time_fil  = psData%i(7);
!  cfl_switch  = psData%i(8);
!  i_time_cont = psData%i(9);
!!----------------------------------------------------------------------------------------
!! terminate if desired
!
!  if(term_status.eq.1) call terminate_run(io,1)
!  return
!!----------------------------------------------------------------------------------------
!! format statements
!
!  1 format(a66,i12)
!  2 format(a66,1pe12.5)
!!----------------------------------------------------------------------------------------
!  return
!  end subroutine read_active_file

!----------------------------------------------------------------------
! New routine by Ramanan. 02/16/05
! This routine will execute the command using either `system' or `ishell'. 
! Gets rid of some of the #if preprocessor commands from the code.
!----------------------------------------------------------------------
  subroutine execute_command(cmd)
  
  character(*), intent(in)  :: cmd

# if SGI || SP2 || CPQ
      call system(trim(cmd))
# endif

# if T3E || X1
      call ishell(trim(cmd))
# endif

  end subroutine execute_command
!----------------------------------------------------------------------


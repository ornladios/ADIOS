!
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!
!  Coupling writer/reader
!
!  Read a 2D array
!  to/from a file or from/to another coupling code using DataSpaces/ADIOS
!
!  nx * ny          processes write a 2D array, where each process writes an
!
!  Data written
!    xy          2D global array as received from the writers
!    xy2         2D array with block,* decomp
!
!    All data are real*8 type
!
! (c) Oak Ridge National Laboratory, 2009
! Author: Norbert Podhorszki
!
module coupling_reader_2D_comm
    ! arguments
    character(len=256) :: filename
    integer :: npx, npy, npz  ! # of processors in x-y-z direction
    integer :: timesteps      ! number of times to read data
    integer :: read_method    ! 0=bp, 1=hdf5, 2=dart, 3=dimes
    integer :: read_mode      ! 0=whole array on each process
                              ! 1=1D decomposition on 1st dim

    integer :: nproc          ! number of processors

    real*8, dimension(:,:),   allocatable :: xy, xy2

    ! Offsets and sizes
    integer :: offs_x, offs_y, offs_z
    integer :: dim_x_local, dim_y_local, dim_z_local
    integer :: dim_x_global, dim_y_global, dim_z_global

    ! MPI variables
    integer :: group_comm, self_comm
    integer :: rank
    integer :: ierr

    ! ADIOS variables
    character (len=200) :: group
    !character (len=6)   :: nprocstr
    integer*8 :: adios_handle, adios_totalsize, adios_groupsize, adios_buf_size
    integer   :: adios_err


    integer   :: ts   ! actual timestep
    integer   :: wts  ! writer's output timestep index (read from 1,2...)
    logical, parameter   :: dump_text = .true.

end module coupling_reader_2D_comm


program coupling
    use coupling_reader_2D_comm
    implicit none
    include 'mpif.h'

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_dup (MPI_COMM_SELF, self_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    call processArgs()
    if (rank == 0) then
        print *, " Output file: "//trim(filename)
        print '(" Process number        : ",i0," x ",i0)',  npx,npy
        print '(" Method for reading: ",i0)',  read_method

        if (nproc .ne. npx*npy) then
            print '(" Error: Number of processors ",i0,"does not match N*M=",i0)', nproc, npx*npy
            call exit(1)
        endif
    endif

    call adios_set_read_method(read_method, ierr) ! 4 = NSSI, 3 = dimes, 2 = dart, 0 = bp
    call adios_read_init (group_comm, ierr)
    call adios_init("coupling2D_reader.xml", adios_err)
    !call MPI_Barrier (group_comm, ierr)

    !write (*,*) 'rank ', rank, "init completed"
    !write (nprocstr,'(i0)') nproc

    call sleep(5)
    call MPI_Barrier (MPI_COMM_WORLD, ierr)
    wts = 1  ! start reading from writer's 1st step
    do ts = 1, timesteps
        call readArrays()
        call printArrays()
        call advanceArrays()
        !print '("rank=",i0," goes to sleep after step ",i0)', rank, ts
        if (ts < timesteps) call sleep(10)
        call MPI_Barrier (MPI_COMM_WORLD, ierr)
        !print '("rank=",i0," woke up")', rank
    enddo

    ! Terminate
    call MPI_Barrier (MPI_COMM_WORLD, ierr)
    call adios_read_finalize (adios_err)
    call adios_finalize (rank, adios_err)
    call MPI_Finalize (ierr)
end program coupling



!!***************************
subroutine readArrays()
    use coupling_reader_2D_comm
    implicit none
    integer             :: i,j,k
    integer             :: gcnt, vcnt, acnt, tstart, ntsteps, tstop, vrank, vtype, time_index
    integer*8           :: fh, gh, read_bytes
    integer*8, dimension(3) :: offset, readsize
    character(len=256)  :: fn

    write (fn,'(a,"_",i3.3,".bp")') trim(filename), wts

    ! Read in data using ADIOS Read API
    call adios_fopen (fh, fn, group_comm, gcnt, adios_err)
    if (adios_err .ne. 0) then
        if (wts .eq. 1) then
            do while (adios_err .ne. 0)
                print '("Writer''s first data ",a," is missing. Wait until it becomes available")', trim(fn)
                call sleep(15)
                call MPI_Barrier (group_comm, ierr)
                call adios_fopen (fh, fn, group_comm, gcnt, adios_err)
            enddo
        else
            print '("Writer data ",a," is missing. Progress without updating reader data")', trim(fn)
            call MPI_Barrier (group_comm, ierr)
            return
        endif
    endif

    wts = wts + 1 ! next time read the next data set

    call adios_gopen (fh, gh, "writer2D", vcnt, acnt, adios_err)
    if (adios_err .ne. 0) then
        call exit(adios_err)
    endif

    offset = 0
    readsize = 1

    ! DEBUG read time_index
    !print '("rank=",i0,": Read in __ADIOS_TIME_INDEX__ from ",a)', rank, trim(fn)
    !call adios_read_var(gh, "__ADIOS_TIME_INDEX__", offset, readsize, time_index, read_bytes)
    !if (read_bytes .lt. 0) then
    !    call exit(read_bytes)
    !elseif (read_bytes .ne. 4) then
    !    print '("Wanted to read __ADIOS_TIME_INDEX__ but read ", i0, " bytes instead of 4")', read_bytes
    !    call exit(read_bytes)
    !endif
    !print '("rank=",i0,": time_index = ",i0)', rank, time_index

    ! read dim_x_global into dim_x_local
    print '("rank=",i0,": Read in dim_x_global from ",a)', rank, trim(fn)
    call adios_read_var(gh, "dim_x_global", offset, readsize, dim_x_global, read_bytes)
    if (read_bytes .lt. 0) then
        call exit(read_bytes)
    elseif (read_bytes .ne. 4) then
        print '("Wanted to read dim_x_global but read ", i0, " bytes instead of 4")', read_bytes
        call exit(read_bytes)
    endif

    ! read dim_y_global into dim_y_local
    print '("rank=",i0,": Read in dim_y_global from ",a)', rank, trim(fn)
    call adios_read_var(gh, "dim_y_global", offset, readsize, dim_y_global, read_bytes)
    if (read_bytes .lt. 0) then
        call exit(read_bytes)
    elseif (read_bytes .ne. 4) then
        print '("Wanted to read dim_y_global but read ", i0, " bytes instead of 4")', read_bytes
        call exit(read_bytes)
    endif

    ! Calculate the local x,y offsets
    if (read_mode == 0) then
        ! each process reads the whole array
        dim_x_local = dim_x_global
        dim_y_local = dim_y_global
        offs_x = 0
        offs_y = 0
    elseif (read_mode == 1) then
        dim_x_local = dim_x_global / nproc
        dim_y_local = dim_y_global
        offs_x = rank * dim_x_local
        offs_y = 0
        ! last process should read all the rest of columns
        if (rank == nproc-1) then
            dim_x_local = dim_x_global - dim_x_local*(nproc-1)
        endif
    endif


    if (.not. allocated(xy)) then
        allocate( xy  (1:dim_x_local, 1:dim_y_local) )
        allocate( xy2 (1:dim_x_local, 1:dim_y_local) )
    endif
    xy = -1.0 ! should be overwritten at reading the array next


    ! Read in the whole xy array
    readsize(1) = dim_x_local
    readsize(2) = dim_y_local
    offset(1)   = offs_x
    offset(2)   = offs_y
    print '("rank=",i0,": Read in xy(",i0,":",i0,",",i0,":",i0,") from ",a)', rank, &
            offset(1)+1, offset(1)+readsize(1), &
            offset(2)+1, offset(2)+readsize(2),trim(fn)
    call adios_read_var(gh, "xy", offset, readsize, xy, read_bytes)
    if (read_bytes .lt. 0) then
        call exit(read_bytes)
    elseif (read_bytes .ne. dim_x_local*dim_y_local*8) then
        print '("Wanted to read xy but read ", i0, " bytes instead of ",i0)', read_bytes, dim_x_local*dim_y_local*8
        call exit(read_bytes)
    endif

    do j=1,dim_y_local
        do i=1,dim_x_local
            xy2(i,j) = rank*1.0 + xy(i,j)/100.0
        enddo
    enddo


    call MPI_Barrier (group_comm, ierr)
    call adios_gclose (gh, adios_err)
    call adios_fclose (fh, adios_err)

end subroutine readArrays

!!***************************
subroutine printArrays()
    use coupling_reader_2D_comm
    implicit none
    integer, parameter  :: u=20
    character(len=256)  :: fn
    integer             :: i,j,k,writer

    if (read_mode == 0) then

        writer = mod(ts,nproc)
        if (writer == rank) then

            write (fn, '("reader_",i3.3,".bp")') ts
            call adios_open (adios_handle, "reader2D", fn, "w", self_comm, adios_err)
#include "gwrite_reader2D.fh"
            ! start streaming from buffer to disk
            call adios_close (adios_handle, adios_err)
            print '("rank=",i0,": Wrote xy to ",a)', rank, trim(fn)

        endif

        if (dump_text) then

            write (fn, '("reader_",i3.3,"_",i0,".txt")') ts, rank
            open (u, FILE=fn, STATUS='NEW', FORM="FORMATTED")

            ! print xyz_bbb
            write (u,'("xy(1:",i0,",1:",i0,") = ")') dim_x_local, dim_y_local
            do j=1,dim_y_local
                do i=1,dim_x_local
                    write (u, '(f6.2," ",$)') xy(i,j)
                enddo
                write (u,*) " "  ! new line
            enddo

            close (u)

        endif

    elseif (read_mode ==1) then
            write (fn, '("reader_",i3.3,".bp")') ts
            call adios_open (adios_handle, "reader2D", fn, "w", group_comm, adios_err)
#include "gwrite_reader2D.fh"
            ! start streaming from buffer to disk
            call adios_close (adios_handle, adios_err)
            print '("rank=",i0,": Wrote xy to ",a)', rank, trim(fn)
    endif

end subroutine printArrays

!!***************************
subroutine advanceArrays()
    use coupling_reader_2D_comm
    implicit none

    xy = 0.9 * xy
    xy2 = 0.9 * xy2

end subroutine advanceArrays


!!***************************
subroutine usage()
    print *, "Usage: coupling_reader_2D file N M method read-mode [timesteps]"
    print *, "file:   name of file to write/read"
    print *, "N:      number of processes in X dimension"
    print *, "M:      number of processes in Y dimension"
    print *, "method: DART or DIMES for memory-to-memory coupling, otherwise file-based"
    print *, "read-mode: 0: each process reads whole global array"
    print *, "           1: 1D decomposition in 1st dimension"
    print *, "timesteps: how many times to write data"
end subroutine usage

!!***************************
subroutine processArgs()
    use coupling_reader_2D_comm

#ifndef __GFORTRAN__
!#ifndef __GNUC__
    interface
         integer function iargc()
         end function iargc
    end interface
!#endif
#endif

    character(len=256) :: npx_str, npy_str, method_str, mode_str, ts_str
    integer :: numargs

    !! process arguments
    numargs = iargc()
    !print *,"Number of arguments:",numargs
    if ( numargs < 5 ) then
        call usage()
        call exit(1)
    endif
    call getarg(1, filename)
    call getarg(2, npx_str)
    call getarg(3, npy_str)
    read (npx_str,'(i5)') npx
    read (npy_str,'(i5)') npy

    call getarg(4, method_str)
    if (trim(method_str) .eq. "DART") then
        read_method = 2
    elseif (trim(method_str) .eq. "DIMES") then
        read_method = 3
    elseif (trim(method_str) .eq. "NSSI") then
        read_method = 4
    else
        read_method = 0
    endif

    call getarg(5, mode_str)
    read (mode_str,'(i5)') read_mode
    if (read_mode < 0 .or. read_mode > 1) then
        print *, "Argument error: read_mode must be 0 or 1"
        call usage()
        call exit(2)
    endif

    if ( numargs == 6 ) then
        call getarg(6, ts_str)
        read (ts_str,'(i6)') timesteps
    else
        timesteps = 1
    endif

end subroutine processArgs

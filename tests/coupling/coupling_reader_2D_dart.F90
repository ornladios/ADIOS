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
!  ndx * ndy        piece with filling with its rank as real*8 value
!  
!  Data written
!    xy         2D array with block,block decomp
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
    integer :: ndx, ndy, ndz  ! size of array per processor
    integer :: timesteps      ! number of times to read data

    integer :: nproc          ! number of processors

    real*8, dimension(:,:),   allocatable :: xy,xy2

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

    ! actual timestep
    integer   :: ts
    logical, parameter   :: dump_text = .true.

    ! DART variables and parameters
    integer, parameter   :: procall = 1024 ! all processes in play (writer+reader+server)


end module coupling_reader_2D_comm


program coupling
    use coupling_reader_2D_comm
    implicit none
    include 'mpif.h'
    integer :: dartrank, dartpeers

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_dup (MPI_COMM_SELF, self_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    print '("rank=",i0," call dart_init(",i0,",2)...")', rank, nproc
    call dart_init (nproc, 2)
    print '("rank=",i0," returned from dart_init()...")', rank
    call adios_init("coupling2D_reader.xml", adios_err)
    !call MPI_Barrier (group_comm, ierr)

    call processArgs()
    if (rank == 0) then
        print *, " Output file: "//trim(filename)
        print '(" Process number        : ",i0," x ",i0)',  npx,npy
        print '(" Array size per process: ",i0," x ",i0)',  ndx,ndy

        if (nproc .ne. npx*npy) then
            print '(" Error: Number of processors ",i0,"does not match N*M=",i0)', nproc, npx*npy
            call exit(1)
        endif
    endif

    !write (*,*) 'rank ', rank, "init completed"
    !write (nprocstr,'(i0)') nproc

    call sleep(5)
    call MPI_Barrier (MPI_COMM_WORLD, ierr)

    print '("rank=",i0," call dart_rank()...")', rank
    call dart_rank(dartrank)
    print '("rank=",i0," call dart_peers()...")', rank
    call dart_peers(dartpeers)
    print '("rank=",i0," returned from dart_peers()...")', rank

    do ts = 1, timesteps
        !call readArraysWithMethod()
        call readArraysWithReadAPI()
        call printArrays()
        !print '("rank=",i0," goes to sleep after step ",i0)', rank, ts
        if (ts < timesteps) call sleep(30)
        call MPI_Barrier (MPI_COMM_WORLD, ierr)
        !print '("rank=",i0," woke up")', rank
    enddo

    ! Terminate
    call MPI_Barrier (MPI_COMM_WORLD, ierr)
    call adios_finalize (rank, adios_err)
    call dart_finalize 
    call MPI_Finalize (ierr)
end program coupling



!!***************************
subroutine readArraysWithReadAPI()
    use coupling_reader_2D_comm
    implicit none
    integer             :: i,j,k
    integer             :: gcnt, vcnt, acnt, tstart, ntsteps, tstop, vrank, vtype, timedim
    integer*8           :: fh, gh, read_bytes
    integer*8, dimension(3) :: offset, readsize
    character(len=256)  :: fn

    call dart_lock_on_read

    ! read dim_x_global into dim_x_local and 
    ! dim_y_global into dim_y_local very first
    if (ts == 1) then
        call dart_get("dim_x_global"//char(0), 0, 4, 0, 0, 0, 0, 0, 0, dim_x_local)
        call dart_get("dim_y_global"//char(0), 0, 4, 0, 0, 0, 0, 0, 0, dim_y_local)
    endif
    
    if (.not. allocated(xy)) then
        allocate( xy(1:dim_x_local, 1:dim_y_local) )
    endif
    xy = -1.0 ! should be overwritten at reading the array next


    ! Read in the whole xy array
    call dart_get("xy"//char(0), 0, 8, 0, 0, 0, dim_y_local-1, dim_x_local-1, 0, xy)
    print '("Get rank=",i0," version=",i0,": xy(1:",i0,",1:",i0,")")', rank, 0, dim_x_local, dim_y_local

    call dart_unlock_on_read

end subroutine readArraysWithReadAPI

!!***************************
subroutine printArrays()
    use coupling_reader_2D_comm
    implicit none
    integer, parameter  :: u=20
    character(len=256)  :: fn
    integer             :: i,j,k,writer

    writer = mod(ts,nproc)
    if (writer == rank) then

        write (fn, '("reader_",i0,".bp")') ts
        call adios_open (adios_handle, "reader2D", fn, "w", self_comm, adios_err)
#include "gwrite_reader2D.fh"
        ! start streaming from buffer to disk
        call adios_close (adios_handle, adios_err)
        print '("rank=",i0,": Wrote xy to ",a)', rank, trim(fn)

    endif

    if (dump_text) then

        write (fn, '("reader_",i0,"_",i0,".txt")') ts, rank
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

end subroutine printArrays

!!***************************
subroutine usage()
    print *, "Usage: genarray file N M nx ny [timesteps]"
    print *, "file:   name of file to write/read"
    print *, "N:      number of processes in X dimension"
    print *, "M:      number of processes in Y dimension"
    print *, "nx:     local array size in X dimension per processor (NOT USED)"
    print *, "ny:     local array size in Y dimension per processor (NOT_USED)"
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

    character(len=256) :: npx_str, npy_str, npz_str, ndx_str, ndy_str, ndz_str, ts_str
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
    !call getarg(4, npz_str)
    read (npx_str,'(i5)') npx
    read (npy_str,'(i5)') npy
    !read (npz_str,'(i5)') npz
    call getarg(4, ndx_str)
    call getarg(5, ndy_str)
    !call getarg(7, ndz_str)
    read (ndx_str,'(i6)') ndx
    read (ndy_str,'(i6)') ndy
    !read (ndz_str,'(i6)') ndz
    if ( numargs == 6 ) then
        call getarg(6, ts_str)
        read (ts_str,'(i6)') timesteps
    else
        timesteps = 1
    endif

end subroutine processArgs

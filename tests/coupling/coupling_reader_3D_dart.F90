!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!
!  Coupling writer/reader
!
!  Write a 3D array 
!  to/from a file or from/to another coupling code using DataSpaces/ADIOS
!
!  nx * ny * nz          processes write a 3D array, where each process writes an
!  ndx * ndy * ndz       piece with filling with its rank as real*8 value
!  
!  Data written
!    xyz         3D array with block,block,block decomp
!
!    All data are real*8 type
!
! (c) Oak Ridge National Laboratory, 2009
! Author: Norbert Podhorszki
!
module coupling_reader_3D_comm
    ! arguments
    character(len=256) :: filename 
    integer :: npx, npy, npz  ! # of processors in x-y-z direction
    integer :: ndx, ndy, ndz  ! size of array per processor
    integer :: timesteps      ! number of times to read data

    integer :: nproc          ! number of processors

    real*8, dimension(:,:,:),   allocatable :: xyz

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


end module coupling_reader_3D_comm


program coupling
    use coupling_reader_3D_comm
    implicit none
    include 'mpif.h'
    integer :: dartrank, dartpeers

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_dup (MPI_COMM_SELF, self_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    print '("rank=",i0," call dart_init()...")', rank
    call dart_init (procall, 2)
    print '("rank=",i0," returned from dart_init()...")', rank
    call adios_init ("coupling3D.xml", ierr)
    !call MPI_Barrier (group_comm, ierr)

    call processArgs()
    if (rank == 0) then
        print *, " Output file: "//trim(filename)
        print '(" Process number        : ",i0," x ",i0," x ",i0)',  npx,npy,npz
        print '(" Array size per process: ",i0," x ",i0," x ",i0)',  ndx,ndy,ndz

        if (nproc .ne. npx*npy*npz) then
            print '(" Error: Number of processors ",i0,"does not match N*M*K=",i0)', nproc, npx*npy*npz
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
    use coupling_reader_3D_comm
    implicit none
    integer             :: i,j,k
    integer             :: gcnt, vcnt, acnt, tstart, ntsteps, tstop, vrank, vtype, timedim
    integer*8           :: fh, gh, read_bytes
    integer*8, dimension(3) :: offset, readsize
    character(len=256)  :: fn

    call dart_lock_on_read

    ! read dim_x_global into dim_x_local and 
    ! dim_y_global into dim_y_local and
    ! dim_z_global into dim_z_local
    if (ts == 1) then
        call dart_get("dim_x_global"//char(0), 0, 4, 0, 0, 0, 0, 0, 0, dim_x_local)
        call dart_get("dim_y_global"//char(0), 0, 4, 0, 0, 0, 0, 0, 0, dim_y_local)
        call dart_get("dim_z_global"//char(0), 0, 4, 0, 0, 0, 0, 0, 0, dim_z_local)
    endif

    if (.not. allocated(xyz)) then
        allocate( xyz(1:dim_x_local, 1:dim_y_local, 1:dim_z_local) )
    endif
    xyz = -1.0 ! should be overwritten at reading the array next


    ! Read in the whole xyz array
    !call dart_get("xyz"//char(0), ts-1, 8, 0, 0, 0, dim_z_local-1, dim_y_local-1, dim_x_local-1, xyz)
    call dart_get("xyz"//char(0), ts-1, 8, 0, 0, 0, dim_y_local-1, dim_x_local-1, dim_z_local-1, xyz)
    print '("Get rank=",i0," version=",i0,": xy(1:",i0,",1:",i0,",1:",i0,")")',  &
           rank, ts-1, dim_x_local, dim_y_local, dim_z_local

    call dart_unlock_on_read

end subroutine readArraysWithReadAPI

!!***************************
subroutine printArrays()
    use coupling_reader_3D_comm
    implicit none
    integer, parameter  :: u=20
    character(len=256)  :: fn
    integer             :: i,j,k,writer

    writer = mod(ts,nproc)
    if (writer == rank) then

        write (fn, '("reader_",i0,".bp")') ts
        call adios_open (adios_handle, "reader3D", fn, "w", self_comm, adios_err)
#include "gwrite_reader3D.fh"
        ! start streaming from buffer to disk
        call adios_close (adios_handle, adios_err)
        print '("rank=",i0,": Wrote xyz to ",a)', rank, trim(fn)

    endif

    if (dump_text) then

        write (fn, '("reader_",i0,"_",i0,".txt")') ts, rank
        open (u, FILE=fn, STATUS='NEW', FORM="FORMATTED")

        ! print xyz
        write (u,'("xyz(1:",i0,",1:",i0,",1:",i0") = ")') dim_x_local, dim_y_local, dim_z_local
        do k=1,dim_z_local
            write (u,'("plane ",i0,":"') k
            do j=1,dim_y_local
                do i=1,dim_x_local
                    write (u, '(f6.2," ",$)') xyz(i,j,k)
                enddo
                write (u,*) " "  ! new line
            enddo
            write (u,*) " "  ! new line
        enddo

        close (u)

    endif

end subroutine printArrays

!!***************************
subroutine usage()
    print *, "Usage: coupling_reader_3D file N M nx ny [timesteps]"
    print *, "file:   name of file to write/read"
    print *, "N:      number of processes in X dimension"
    print *, "M:      number of processes in Y dimension"
    print *, "nx:     local array size in X dimension per processor"
    print *, "ny:     local array size in Y dimension per processor"
    print *, "timesteps: how many times to write data"
end subroutine usage

!!***************************
subroutine processArgs()
    use coupling_reader_3D_comm

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
    if ( numargs < 7 ) then
        call usage()
        call exit(1)
    endif
    call getarg(1, filename)
    call getarg(2, npx_str)
    call getarg(3, npy_str)
    call getarg(4, npz_str)
    read (npx_str,'(i5)') npx
    read (npy_str,'(i5)') npy
    read (npz_str,'(i5)') npz
    call getarg(5, ndx_str)
    call getarg(6, ndy_str)
    call getarg(7, ndz_str)
    read (ndx_str,'(i6)') ndx
    read (ndy_str,'(i6)') ndy
    read (ndz_str,'(i6)') ndz
    if ( numargs == 8 ) then
        call getarg(8, ts_str)
        read (ts_str,'(i6)') timesteps
    else
        timesteps = 1
    endif

end subroutine processArgs

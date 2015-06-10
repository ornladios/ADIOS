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
module coupling_writer_3D_comm
    ! arguments
    character(len=256) :: filename 
    integer :: npx, npy, npz  ! # of processors in x-y-z direction
    integer :: ndx, ndy, ndz  ! size of array per processor
    integer :: timesteps      ! number of timesteps to run

    integer :: nproc          ! number of processors

    real*8, dimension(:,:,:),   allocatable :: xyz

    ! Offsets and sizes
    integer :: offs_x, offs_y, offs_z
    integer :: dim_x_local, dim_y_local, dim_z_local
    integer :: dim_x_global, dim_y_global, dim_z_global

    ! MPI variables
    integer :: group_comm
    integer :: rank
    integer :: ierr

    ! ADIOS variables
    character (len=200) :: group
    !character (len=6)   :: nprocstr
    integer*8 :: adios_handle, adios_totalsize, adios_groupsize, adios_buf_size
    integer   :: adios_err

    ! actual timestep
    integer   :: ts 

    ! DART variables and parameters
    integer, parameter   :: procall = 1024 ! all processes in play (writer+reader+server)

end module coupling_writer_3D_comm


program coupling
    use coupling_writer_3D_comm
    implicit none
    include 'mpif.h'
    integer :: t
    integer :: dartrank, dartpeers

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    print '("rank=",i0," call dart_init()...")', rank
    call dart_init (procall, 1)
    print '("rank=",i0," returned from dart_init()")', rank
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

    print '("rank=",i0," call dart_rank()...")', rank
    call dart_rank(dartrank)
    print '("rank=",i0," call dart_peers()...")', rank
    call dart_peers(dartpeers)
    print '("rank=",i0," returned from dart_peers()...")', rank

    ! Calculate global size
    call allocateLocalArrays()
    do ts=1,timesteps
        call generateLocalArrays()
        call writeArrays()
        call MPI_Barrier (MPI_COMM_WORLD, ierr)
        !print '("rank=",i0," goes to sleep after step ",i0)', rank, ts
        if (ts < timesteps) call sleep(30)
        !print '("rank=",i0," woke up")', rank
    enddo

    ! Terminate
    call MPI_Barrier (MPI_COMM_WORLD, ierr)
    call dart_finalize 
    call MPI_Finalize (ierr)
end program coupling


!!***************************
subroutine allocateLocalArrays()
    use coupling_writer_3D_comm
    implicit none
    integer :: posx, posy, posz ! position index in the array
    integer :: i,j,k


!!  3D array with block,block decomposition
    posx = mod(rank, npx)     ! 1st dim easy: 0, npx, 2npx... are in the same X position
    posy = mod(rank/npx, npy) ! 2nd dim: (0, npx-1) have the same dim (so divide with npx first)
    posz = rank/(npx*npy)     ! 3rd dim: npx*npy processes belong into one dim
    offs_x = posx * ndx
    offs_y = posy * ndy
    offs_z = posz * ndz
    dim_x_local   = ndx
    dim_y_local   = ndy
    dim_z_local   = ndz
    dim_x_global  = npx * ndx
    dim_y_global  = npy * ndy
    dim_z_global  = npz * ndz

    print '("rank=",i0,": 3D array pos: ",i0,",",i0,",",i0," offset: ",i0,",",i0,",",i0)',  &
             rank, posx, posy, posz, offs_x, offs_y, offs_z

    allocate( xyz(1:ndx, 1:ndy, 1:ndz) )

end subroutine allocateLocalArrays

!!***************************
subroutine generateLocalArrays()
    use coupling_writer_3D_comm
    implicit none

    print '("rank=",i0," set matrix to ",f6.2)', rank, 1.0*rank+0.01*ts
    xyz = 1.0*rank + 0.01*ts

end subroutine generateLocalArrays

!!***************************
subroutine writeArrays()
    use coupling_writer_3D_comm
    implicit none
    character(len=256) :: fn
    ! Write out data using DART
    ! ??? Why do we need this ???
    call dart_lock_on_write
    write (fn,'(a,"_",i0,".bp")') trim(filename), ts
    !print '("rank=",i0," file=",A )', rank,  trim(filename)

    !! Syntax is: var_name, version, size_elem, {bounding box}, data values.
    !! bounding box: {slowest changing dimension, ..., fastest changing dimension}
    !!               { start, end,  ..., start, end}
    !! call dart_put("m3d"//char(0), i-1, 8, 0, j-1, 0, n-1, k-1, n-1, m3d(j:k,1:n,1:n))
    !! call dart_put("m3d"//char(0), i-1, 8, 0, j-1, 0, n-1, k-1, n-1, m3d)    
    if (rank == 0 .and. ts == 1) then
        print '("Put dim_x_global = ",i0," dim_y_global = ",i0," dim_z_global = ",i0," version 0")', &
                dim_x_global, dim_y_global, dim_z_global
        call dart_put("dim_x_global"//char(0), 0, 4, 0, 0, 0, 0, 0, 0, dim_x_global)
        call dart_put("dim_y_global"//char(0), 0, 4, 0, 0, 0, 0, 0, 0, dim_y_global)
        call dart_put("dim_z_global"//char(0), 0, 4, 0, 0, 0, 0, 0, 0, dim_z_global)
    endif
    print '("Put rank=",i0," version=",i0,": xyz(",i0,":",i0,",",i0,":",i0,",",i0,":",i0")")', &
            rank, ts-1, offs_x, offs_x+dim_x_local, offs_y, offs_y+dim_y_local, offs_z, offs_z+dim_z_local
    !call dart_put("xyz"//char(0), ts-1, 8, &
    !              offs_z, offs_y, offs_x, &
    !              offs_z+dim_z_local-1, offs_y+dim_y_local-1, offs_x+dim_x_local-1, &
    !              xyz)    
    call dart_put("xyz"//char(0), ts-1, 8, &
                  offs_y, offs_x, offs_z, &
                  offs_y+dim_y_local-1, offs_x+dim_x_local-1, offs_z+dim_z_local-1, &
                  xyz)    

    ! ??? Why do we need this ???
    call dart_put_sync
    
    print '("rank=",i0,": write completed")', rank
    call dart_unlock_on_write
end subroutine writeArrays


!!***************************
subroutine usage()
    print *, "Usage: coupling_writer_3D file N M K nx ny nz [timesteps]"
    print *, "file:   name of file to write/read"
    print *, "N:      number of processes in X dimension"
    print *, "M:      number of processes in Y dimension"
    print *, "K:      number of processes in Z dimension"
    print *, "nx:     local array size in X dimension per processor"
    print *, "ny:     local array size in Y dimension per processor"
    print *, "nz:     local array size in Z dimension per processor"
    print *, "timesteps: how many times to write data"
end subroutine usage

!!***************************
subroutine processArgs()
    use coupling_writer_3D_comm

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

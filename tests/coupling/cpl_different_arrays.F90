!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!
!  Coupling writer/reader
!
!  Write/read different kind of arrays and decompositions in several steps
!  to/from a file or from/to another coupling code using DataSpaces/ADIOS
!
!  nx * ny * nz     processes write a 3D array, where each process writes an
!  ndx * ndy * ndz  piece with filling with its rank as real*8 value
!  
!  Data written
!    xyz_bbb    3D array with 3D (block,block,block) decomposition
!    xyz_bss    3D array with 1D (block,*,*) decomposition
!    !  xyz_css    3D array with 1D cyclic (cyclic,*,*) decomposition
!    xy_bb      2D array with block,block decomp
!    xy_bs      2D array with block,* decomp
!    !  xy_cs      2D array with cyclic,* decomp
!    x_b        1D array with block decomp
!    !  x_c        1D array with cyclic decomp
!    scalar     A scalar value
!
!    All data are real*8 type
!
! (c) Oak Ridge National Laboratory, 2009
! Author: Norbert Podhorszki
!
module coupling_comm
    ! arguments
    character(len=256) :: filename 
    integer :: npx, npy, npz  ! # of processors in x-y-z direction
    integer :: ndx, ndy, ndz  ! size of array per processor
    logical :: do_write       ! .true.  write the arrays
                              ! .false. read the arrays and then write them to file

    integer :: nproc          ! number of processors

    !real*8, dimension(:,:,:), allocatable :: xyz
    real*8, dimension(:,:,:), allocatable :: xyz_bbb, xyz_bss
    real*8, dimension(:,:),   allocatable :: xy_bb, xy_bs
    real*8, dimension(:),     allocatable :: x_b

    ! Offsets and sizes
    integer :: offs_x_b, dim_x_b                                ! x_b
    integer :: offs_x_bb, offs_y_bb, dim_x_bb, dim_y_bb         ! xy_bb
    integer :: offs_x_bs, dim_x_bs                              ! xy_bs
    integer :: offs_x_bbb, offs_y_bbb, offs_z_bbb               ! xyz_bbb
    integer :: dim_x_bbb, dim_y_bbb, dim_z_bbb
    integer :: offs_x_bss, dim_x_bss                            ! xyz_bss

    ! MPI variables
    integer :: group_comm
    integer :: rank
    integer :: ierr

    ! ADIOS variables
    character (len=200) :: group
    !character (len=6)   :: nprocstr
    integer*8 :: adios_handle, adios_totalsize, adios_groupsize, adios_buf_size
    integer   :: adios_err

end module coupling_comm


program coupling
    use coupling_comm
    implicit none
    include 'mpif.h'

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    call adios_init ("cpl_different_arrays.xml", ierr)
    !call MPI_Barrier (group_comm, ierr)

    call processArgs()
    if (rank == 0) then
        if (do_write) then
            print *, " Output file: "//trim(filename)
        else
            print *, " Input file: "//trim(filename)
        endif
        print '(" Process number        : ",i0," x ",i0," x ",i0)',  npx,npy,npz
        print '(" Array size per process: ",i0," x ",i0," x ",i0)',  ndx,ndy,ndz

        if (nproc .ne. npx*npy*npz) then
            print '(" Error: Number of processors ",i0,"does not match ndx*ndy*ndz=",i0)', nproc, npx*npy*npz
            call exit(1)
        endif
    endif

    !write (*,*) 'rank ', rank, "init completed"
    !write (nprocstr,'(i0)') nproc

    ! Calculate global size
    call generateLocalArrays()
    if (do_write) then
        call writeArrays()
    else
        call readArrays()
        call printArrays()
    endif

    ! Terminate
    call MPI_Barrier (MPI_COMM_WORLD, ierr)
    call adios_finalize (rank, adios_err)
    call MPI_Finalize (ierr)
end program coupling


!!***************************
subroutine generateLocalArrays()
    use coupling_comm
    implicit none
    integer :: posx, posy, posz ! position index in the array
    integer :: i,j,k

!!  3D array with block,block,block decomposition
    posx = mod(rank, npx)     ! 1st dim easy: 0, npx, 2npx... are in the same X position
    posy = mod(rank/npx, npy) ! 2nd dim: (0, npx-1) have the same dim (so divide with npx first)
    posz = rank/(npx*npy)     ! 3rd dim: npx*npy processes belong into one dim
    offs_x_bbb = posx * ndx
    offs_y_bbb = posy * ndy
    offs_z_bbb = posz * ndz
    dim_x_bbb  = npx * ndx
    dim_y_bbb  = npy * ndy
    dim_z_bbb  = npz * ndz

    print '("3D bbb rank=",i0," pos: ",i0,",",i0,",",i0," offset: ",i0,",",i0,","i0)',  &
             rank, posx, posy, posz, offs_x_bbb, offs_y_bbb, offs_z_bbb
    allocate( xyz_bbb(1:ndx, 1:ndy, 1:ndz) )
    do k=1,ndz
        do j=1,ndy
            do i=1,ndx
                xyz_bbb(i,j,k) = 1.0*rank 
            enddo
        enddo
    enddo

!!  3D array with block,*,* decomposition
    offs_x_bss = rank * ndx
    dim_x_bss  = nproc * ndx
    print '("3D bss rank=",i0," offset: ",i0," 0,0")', rank, offs_x_bss
    allocate( xyz_bss(1:ndx, 1:ndy, 1:ndz) )
    do k=1,ndz
        do j=1,ndy
            do i=1,ndx
                xyz_bss(i,j,k) = 1.0*rank 
            enddo
        enddo
    enddo


!!  2D array with block,block decomposition
    posx = mod(rank, npx)     ! 1st dim easy: 0, npx, 2npx... are in the same X position
    posy = rank/npx           ! 2nd dim: npx processes belong into one dim
    offs_x_bb = posx * ndx
    offs_y_bb = posy * ndy
    dim_x_bb  = npx * ndx
    dim_y_bb  = npy * npz * ndy

    print '("2D bb rank=",i0," pos: ",i0,",",i0," offset: ",i0,",",i0)',  &
             rank, posx, posy, offs_x_bb, offs_y_bb

    allocate( xy_bb(1:ndx, 1:ndy) )
    do j=1,ndy
        do i=1,ndx
            xy_bb(i,j) = 1.0*rank 
        enddo
    enddo


!!  2D array with block,*,* decomposition
    offs_x_bs = rank * ndx
    dim_x_bs  = nproc * ndx
    print '("2D bs rank=",i0," offset: ",i0," 0,0")', rank, offs_x_bs
    allocate( xy_bs(1:ndx, 1:ndy) )
    do j=1,ndy
        do i=1,ndx
            xy_bs(i,j) = 1.0*rank 
        enddo
    enddo
 
!!  1D array with block decomposition
    offs_x_b = rank * ndx
    dim_x_b  = nproc * ndx
    print '("1D bs rank=",i0," offset: ",i0," 0,0")', rank, offs_x_b
    allocate( x_b(1:ndx) )
    x_b = 1.0*rank

end subroutine generateLocalArrays


!!***************************
subroutine writeArrays()
    use coupling_comm
    implicit none
    ! Write out data using ADIOS

    !print '("rank=",i0," file=",A )', rank,  trim(filename)
    call adios_open (adios_handle, "coupling", filename, "w", group_comm, adios_err)

    adios_groupsize =   24 * 4                & ! bunch of integers
                + 8 * (ndx) * (ndy)           & ! xy_bb
                + 8 * (ndx) * (ndy) * (ndz)     ! xyz_bbb
               ! + 8 * (ndx) &                   ! x_b
               ! + 8 * (ndx) * (ndy) &           ! xy_bs
               ! + 8 * (ndx) * (ndy) * (ndz)     ! xyz_bss

#include "gwrite_coupling.fh"

    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)
    call adios_write (adios_handle, "ndx", ndx, adios_err)
    call adios_write (adios_handle, "ndy", ndy, adios_err)
    call adios_write (adios_handle, "ndz", ndz, adios_err)
    call adios_write (adios_handle, "nproc", nproc, adios_err)
    call adios_write (adios_handle, "npx", npx, adios_err)
    call adios_write (adios_handle, "npy", npy, adios_err)
    call adios_write (adios_handle, "npz", npz, adios_err)
    !call adios_write (adios_handle, "offs_x_b", offs_x_b, adios_err)
    !call adios_write (adios_handle, "dim_x_b", dim_x_b, adios_err)
    call adios_write (adios_handle, "offs_x_bb", offs_x_bb, adios_err)
    call adios_write (adios_handle, "offs_y_bb", offs_y_bb, adios_err)
    call adios_write (adios_handle, "dim_x_bb", dim_x_bb, adios_err)
    call adios_write (adios_handle, "dim_y_bb", dim_y_bb, adios_err)
    !call adios_write (adios_handle, "offs_x_bs", offs_x_bs, adios_err)
    !call adios_write (adios_handle, "dim_x_bs", dim_x_bs, adios_err)
    call adios_write (adios_handle, "offs_x_bbb", offs_x_bbb, adios_err)
    call adios_write (adios_handle, "offs_y_bbb", offs_y_bbb, adios_err)
    call adios_write (adios_handle, "offs_z_bbb", offs_z_bbb, adios_err)
    call adios_write (adios_handle, "dim_x_bbb", dim_x_bbb, adios_err)
    call adios_write (adios_handle, "dim_y_bbb", dim_y_bbb, adios_err)
    call adios_write (adios_handle, "dim_z_bbb", dim_z_bbb, adios_err)
    !call adios_write (adios_handle, "offs_x_bss", offs_x_bss, adios_err)
    !call adios_write (adios_handle, "dim_x_bss", dim_x_bss, adios_err)
    !call adios_write (adios_handle, "x_b", x_b, adios_err)
    call adios_write (adios_handle, "xy_bb", xy_bb, adios_err)
    !call adios_write (adios_handle, "xy_bs", xy_bs, adios_err)
    call adios_write (adios_handle, "xyz_bbb", xyz_bbb, adios_err)
    !call adios_write (adios_handle, "xyz_bss", xyz_bss, adios_err)

    ! start streaming from buffer to disk
    call adios_close (adios_handle, adios_err)
    print '("rank=",i0,": write completed")', rank
end subroutine writeArrays


!!***************************
subroutine readArrays()
    use coupling_comm
    implicit none
    integer, parameter  :: u=10
    character(len=256)  :: fn
    integer             :: i,j,k
    ! Read in data using ADIOS
    call adios_open (adios_handle, "coupling", filename, "r", group_comm, adios_err)

#include "gread_coupling.fh"

    ! start streaming from disk to buffer 
    call adios_close (adios_handle, adios_err)
    
    write (fn, '("data",i0,".txt")') rank
    open (u, FILE=fn, STATUS='NEW', FORM="FORMATTED")

    ! print xyz_bbb
    write (u,'("xyz_bbb(",i0,":",i0,":",i0,") = ")') ndx, ndy, ndz
    do k=1,ndz
        do j=1,ndy
            do i=1,ndx
                write (u, '(f6.2," ",$)') xyz_bbb(i,j,k)
            enddo
            write (u,*) " "  ! new line
        enddo
    enddo

    close (u)

end subroutine readArrays


!!***************************
subroutine printArrays()
    use coupling_comm
    implicit none
end subroutine printArrays


!!***************************
subroutine usage()
    print *, "Usage: genarray mode file N  M  K  nx  ny  nz"
    print *, "mode:   w | r   to write or read "
    print *, "file:   name of file to write/read"
    print *, "N:      number of processes in X dimension"
    print *, "M:      number of processes in Y dimension"
    print *, "K:      number of processes in Z dimension"
    print *, "nx:     local array size in X dimension per processor"
    print *, "ny:     local array size in Y dimension per processor"
    print *, "nz:     local array size in Z dimension per processor"
end subroutine usage

!!***************************
subroutine processArgs()
    use coupling_comm

#ifndef __GFORTRAN__
!#ifndef __GNUC__
    interface
         integer function iargc()
         end function iargc
    end interface
!#endif
#endif

    character(len=256) :: npx_str, npy_str, npz_str, ndx_str, ndy_str, ndz_str 
    character(len=256) :: mode_str
    integer :: numargs

    !! process arguments
    numargs = iargc()
    !print *,"Number of arguments:",numargs
    if ( numargs < 5 ) then
        call usage()
        call exit(1)
    endif
    call getarg(1, mode_str)
    call getarg(2, filename)
    call getarg(3, npx_str)
    call getarg(4, npy_str)
    call getarg(5, npz_str)
    read (npx_str,'(i5)') npx
    read (npy_str,'(i5)') npy
    read (npz_str,'(i5)') npz
    if (mode_str == "w") then
        do_write = .true.
    else if (mode_str == "r") then 
        do_write = .false.
    else
        print *, "First argument must be w or r"
        call usage()
        call exit(1)
    endif
    call getarg(6, ndx_str)
    call getarg(7, ndy_str)
    call getarg(8, ndz_str)
    read (ndx_str,'(i6)') ndx
    read (ndy_str,'(i6)') ndy
    read (ndz_str,'(i6)') ndz

end subroutine processArgs

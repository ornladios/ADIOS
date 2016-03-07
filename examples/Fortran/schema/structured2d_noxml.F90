! 
 ! ADIOS is freely available under the terms of the BSD license described
 ! in the COPYING file in the top level directory of this source distribution.
 !
 ! Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

! ADIOS C Example: write variables along with an unstructured mesh. 

subroutine printUsage (prgname)
    implicit none
    character(*),   intent(in)  :: prgname
    print *, "Usage: mpirun -np <N> ", prgname, "<nx> <ny>"
    print *, "    <nx> <ny>  2D decomposition values in each dimension of an 2D array"
    print *, "         The product of these number must be equal the number of processes"
    print *, "         e.g. for N=12 you may use  4 3"
end subroutine

subroutine processArgs
!    integer,        intent(out)  :: nproc, npx, npy
    integer                      :: nn
    character(len=32)            :: arg
    COMMON nproc, npx, npy, arg_err

    nn = command_argument_count()
!    write (*,*) 'number of command arguments = ', nn
!    write (*,*) 'number of processor = ', nproc
    if (nn /= 2) then
        call getarg(0, arg)
        call printUsage (arg)
        arg_err = 1
     else
        call getarg(1, arg)  ! write (arg, '(i10)') nproc    
        read (arg, *) npx
!        write (*,*) 'number of x processor = ', npx
        call getarg(2, arg)
        read (arg, *) npy
!        write (*,*) 'number of y processor = ', npy
     endif

     if (npx*npy /= nproc) then
        arg_err = 1
        print *, "ERROR: Product of decomposition numbers in X and Y dimension ", npx*npx
        print *, " != number of processes ", nproc
     endif
end subroutine

program structured2d_f_noxml
    use adios_write_mod
    implicit none
    include 'mpif.h'
    character(len=256)      :: filename = "structured2d_f_noxml.bp"
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle
    integer*8               :: m_adios_group
    integer*8               :: varid
    integer*4               :: comm
    real*8, dimension(:), allocatable       :: X   ! X coordinate
    real*8, dimension(:), allocatable       :: Y   ! Y coordinate 
    real*8, dimension(:), allocatable       :: data
    character(len=20)                       :: schema_version, dimemsions

    integer*4               :: rank, i, j, ierr, adios_err
    integer*4               :: ndx, ndy                            ! size of array per processor
    integer*4               :: offs_x, offs_y                      ! offset in x and y direction
    integer*4               :: nx_local, ny_local                  ! local address
    integer*4               :: posx, posy                          ! position index in the array
    integer*4               :: nx_global, ny_global                ! global address

    !will work with 12 cores, which are arranged by npx=4, npy=3 (4x3)
    integer*4               :: npx                                 ! # of procs in x direction
    integer*4               :: npy                                 ! # of procs in y direction
    integer*4               :: nproc                               ! # of total procs
    integer*4               :: arg_err

    COMMON nproc, npx, npy, arg_err

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, nproc, ierr)

    arg_err = 0
    call processArgs

    if (arg_err .EQ. 1) then
        stop
    endif

    ! will work with 12 cores, which are arranged by npx=4, npy=3 (4x3)
    ndx = 65
    ndy = 129

    !2D array with block,block decomposition
    posx = mod(rank,npx)           ! 1st dim
    posy = rank/npx;               ! 2nd dim
    offs_x = posx * ndx
    offs_y = posy * ndy
    nx_local = ndx
    ny_local = ndy
    nx_global = npx * ndx
    ny_global = npy * ndy

    schema_version = "1.1"
    dimemsions = "nx_global,ny_global"

    allocate (data(0:(ndx-1)*(ndy-1)*8))
    do i = 0, ndx-1
        do j = 0, ndy-1
            data(i*ndy+j) = 1.0*DBLE(rank)
        enddo
    enddo

    allocate (X(0:(ndx-1)*(ndy-1)*8))
    do i = 0, ndx-1
        do j = 0, ndy-1
            X(i*ndy+j) = offs_x + posy*ndx + i*ndx/ndx + DBLE(ndx)*j/DBLE(ndy)
        enddo
    enddo

    allocate (Y(0:(ndx-1)*(ndy-1)*8))
    Y(0) = offs_y
    do i = 0, ndx-1
        do j = 0, ndy-1
            Y(i*ndy+j) = offs_y + ndy*j/ndy 
        enddo
    enddo

    call adios_init_noxml (comm, adios_err)
    call adios_set_max_buffer_size (50) 
    call adios_declare_group (m_adios_group, "structured2d", "", 1, adios_err)
    call adios_select_method (m_adios_group, "MPI", "", "", adios_err)

    ! This example doesn't use varid during writing.
    ! So we simply put 'varid' everywhere.
    ! define a integer
    call adios_define_var (m_adios_group, "nx_global" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (m_adios_group, "ny_global" &
                ,"", adios_integer &
                ,"" ,"","", varid)
    call adios_define_var (m_adios_group, "nproc" &
                ,"", adios_integer &
                ,"" ,"","", varid)
    call adios_define_var (m_adios_group, "offs_x" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (m_adios_group, "offs_y" &
                ,"", adios_integer &
                ,"","" ,"" , varid)
    call adios_define_var (m_adios_group, "nx_local" &
                ,"", adios_integer &
                ,"","" , "", varid)
    call adios_define_var (m_adios_group, "ny_local" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (m_adios_group, "X" &
                    ,"", adios_double &
                    ,"ny_local,nx_local", "ny_global,nx_global", "offs_y,offs_x", varid)
    call adios_define_var (m_adios_group, "Y" &
                    ,"", adios_double &
                    ,"ny_local,nx_local", "ny_global,nx_global", "offs_y,offs_x", varid)
    call adios_define_var (m_adios_group, "data" &
                    ,"", adios_double &
                    ,"ny_local,nx_local", "ny_global,nx_global", "offs_y,offs_x", varid)

    call adios_define_schema_version (m_adios_group, schema_version)
    call adios_define_mesh_structured (dimemsions, "X,Y", "2", m_adios_group, "structuredmesh")
    call adios_define_mesh_timevarying ("no", m_adios_group, "structuredmesh")
    call adios_define_var_mesh (m_adios_group, "data", "structuredmesh")
    call adios_define_var_centering (m_adios_group, "data", "point")

    adios_groupsize = 7*4 & !int
    + 3 * 8 * (nx_local*ny_local) ! double (data)

    call adios_open (adios_handle, "structured2d", filename, "w", comm, adios_err)
    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

    call adios_write (adios_handle, "nx_global", nx_global, adios_err)
    call adios_write (adios_handle, "ny_global", ny_global, adios_err)
    call adios_write (adios_handle, "nproc", nproc, adios_err)
    call adios_write (adios_handle, "offs_x", offs_x, adios_err)
    call adios_write (adios_handle, "offs_y", offs_y, adios_err)
    call adios_write (adios_handle, "nx_local", nx_local, adios_err)
    call adios_write (adios_handle, "ny_local", ny_local, adios_err)

    call adios_write (adios_handle, "X", X, adios_err)
    call adios_write (adios_handle, "Y", Y, adios_err)
    call adios_write (adios_handle, "data", data, adios_err)

    call adios_close (adios_handle, adios_err)
    call MPI_Barrier (comm, ierr)
    call adios_finalize (rank, adios_err)
    call MPI_Finalize (ierr)

end program


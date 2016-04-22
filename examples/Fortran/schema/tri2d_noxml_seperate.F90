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

program tri2d_f_noxml
    use adios_write_mod
    implicit none
    include 'mpif.h'
!    character(len=256)      :: filename = "tri2d_f_noxml.bp"
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: data_file
    integer*8               :: mesh_file
    integer*8               :: data_group
    integer*8               :: mesh_group
    integer*8               :: varid
    integer*4               :: comm

    integer*4               :: rank, i, j, p, p1, ierr, adios_err
    integer*4               :: npoints, num_cells
    integer*4               :: ndx, ndy                            ! size of array per processor
    real*8, dimension(:), allocatable       :: N                                   ! node centered variable
    real*8, dimension(:), allocatable       :: C                                   ! cell centered variable
    real*8, dimension(:,:), allocatable     :: points                              ! X,Y coordinate
    integer*4, dimension(:,:), allocatable  :: cells
    character(len=20)                       :: schema_version, dimemsions

    integer*4               :: offs_x, offs_y                      ! offset in x and y direction
    integer*4               :: nx_local, ny_local                  ! local address
    integer*4               :: posx, posy                          ! position index in the array
    integer*4               :: nx_global, ny_global                ! global address
    integer*4               :: lp, op, lc, oc
    
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

    if (arg_err == 1) then
        stop
    endif

    ndx = 4
    ndy = 3

    npoints = ndx * ndy * npx * npy
    num_cells = (ndx * npx - 1) * (ndy * npy -1) * 2

! 2D array with block,block decomposition
    posx = mod(rank,npx)           ! 1st dim
    posy = rank/npx;               ! 2nd dim
    offs_x = posx * ndx
    offs_y = posy * ndy
    nx_local = ndx
    ny_local = ndy
    nx_global = npx * ndx
    ny_global = npy * ndy

! local mesh + data
    allocate (N(0:ndx*ndy*8-1))            ! 8 is sizeof(double)
    allocate (points(2,0:ndx*ndy*8-1))   ! 8 is sizeof(double)

    if (posx == npx-1 .AND. posy < npy-1 ) then
        allocate (C(0:8*(ndx-1)*ndy*2-1))      ! 8 is sizeof(double)
        allocate (cells(3,0:4*(ndx-1)*ndy*2-1))! 4 is sizeof(int)
    else if (posy == npy-1 .AND. posx < npx-1 ) then
        allocate (C(0:8*ndx*(ndy-1)*2-1))
        allocate (cells(3,0:4*ndx*(ndy-1)*2-1))
    else if ( posx == npx-1 .AND. posy == npy-1 ) then 
        allocate (C(0:8*(ndx-1)*(ndy-1)*2-1))
        allocate (cells(3,0:4*(ndx-1)*(ndy-1)*2-1))
    else
        allocate (C(0:8*ndx*ndy*2-1))
        allocate (cells(3,0:4*ndx*ndy*2-1))
    endif

! generate local data
    lp = ndx * ndy                         ! number of points in this local array
    op = ndx * ndy * rank                  ! offset in global points array

    do i = 0, ndx-1
        do j = 0, ndy-1
            points(1,(i*ndy + j)) = offs_x + posy*ndx + DBLE(i)*ndx/ndx + DBLE(ndx)*j/ndy
            points(2,(i*ndy + j)) = offs_y + DBLE(ndy)*j/ndy
        enddo
    enddo

    do i = 0, lp-1
        N(i) = 1.0*DBLE(rank)
    enddo
    
    if( posx == npx-1 .AND. posy < npy-1 ) then
        lc = (ndx-1)*ndy*2
        oc = posx*ndx*ndy*2 + posy*(ndx*npx-1)*ndy*2
        !write (*,'("Rank",i0,": case 1: define ",i0," triangles, oc=",i0)'), rank, lc, oc
        do i = 0, ndx-2
            do j = 0, ndy-1
               p = i*ndy+j
               if( i<ndx-1 .AND. j<ndy-1 ) then
                   cells(:,2*p)   = (/ op+p, op+p+ndy,   op+p+ndy+1  /)
                   cells(:,2*p+1) = (/ op+p, op+p+ndy+1, op+p+1      /)
                   !write (*,'("Rank",i0,": 1A cell ",i0,"= [",3i3,"]")'), rank, 2*p, cells(:,2*p)
                   !write (*,'("Rank",i0,": 1A cell ",i0,"= [",3i3,"]")'), rank, 2*p+1, cells(:,2*p+1)
                else                     ! extend in Y direction only
                   cells(:,2*p)   = (/ op+p, op+p+ndy,   op+nx_global*ndy+(i+1)*ndy  /)
                   cells(:,2*p+1) = (/ op+p, op+nx_global*ndy+(i+1)*ndy, op+nx_global*ndy+i*ndy /)
                   !write (*,'("Rank",i0,": 1B cell ",i0,"= [",3i3,"]")'), rank, 2*p, cells(:,2*p)
                   !write (*,'("Rank",i0,": 1B cell ",i0,"= [",3i3,"]")'), rank, 2*p+1, cells(:,2*p+1)
                endif
            enddo
        enddo
    else if( posy == npy-1 .AND. posx < npx-1 ) then
        lc = ndx*(ndy-1)*2
        oc = posy*(ndx*npx-1)*ndy*2 + posx*ndx*(ndy-1)*2
        !write (*,'("Rank",i0,": case 2: define ",i0," triangles, oc=",i0)'), rank, lc, oc
        do i = 0, ndx-1
            do j = 0, ndy-2
                p = i*(ndy-1)+j
                p1 = i*ndy+j
                if( i<ndx-1 .AND. j<ndy-1 ) then
                   cells(:,2*p)   = (/ op+p1, op+p1+ndy,   op+p1+ndy+1  /)
                   cells(:,2*p+1) = (/ op+p1, op+p1+ndy+1, op+p1+1      /)
                else
                   cells(:,2*p)   = (/ op+p1, op+ndx*ndy+j,   op+ndx*ndy+j+1  /)
                   cells(:,2*p+1) = (/ op+p1, op+ndx*ndy+j+1, op+p1+1      /)
                endif
            enddo
        enddo
    else if( posx == npx-1 .AND. posy == npy-1 ) then
        lc = (ndx-1)*(ndy-1)*2
        oc = posy*(ndx*npx-1)*ndy*2 + posx*ndx*(ndy-1)*2
        !write (*,'("Rank",i0,": case 3: define ",i0," triangles, oc=",i0)'), rank, lc, oc
        do i = 0, ndx-2
            do j = 0, ndy-2
                p = i*(ndy-1)+j
                p1 = i*ndy+j
                cells(:,2*p)   = (/ op+p1, op+p1+ndy,   op+p1+ndy+1  /)
                cells(:,2*p+1) = (/ op+p1, op+p1+ndy+1, op+p1+1      /)
            enddo
        enddo
    else
        lc = ndx*ndy*2
        oc = posx*ndx*ndy*2 + posy*(ndx*npx-1)*ndy*2
        !write (*,'("Rank",i0,": case 4: define ",i0," triangles, oc=", i0)'), rank, lc, oc
        do i = 0, ndx-1
            do j = 0, ndy-1
                p = i*ndy+j
                if( i<ndx-1 .AND. j<ndy-1 ) then
                    cells(:,2*p)   = (/ op+p, op+p+ndy,   op+p+ndy+1  /)
                    cells(:,2*p+1) = (/ op+p, op+p+ndy+1, op+p+1      /)
                else if( i==ndx-1 .AND. j<ndy-1 ) then
                    cells(:,2*p)   = (/ op+p, op+ndx*ndy+j,   op+ndx*ndy+j+1  /)
                    cells(:,2*p+1) = (/ op+p, op+ndx*ndy+j+1, op+p+1  /)
                else if( i<ndx-1 .AND. j==ndy-1 ) then
                    cells(:,2*p)   = (/ op+p, op+p+ndy,   op+nx_global*ndy+(i+1)*ndy /)
                    cells(:,2*p+1) = (/ op+p, op+nx_global*ndy+(i+1)*ndy, op+nx_global*ndy+i*ndy /)
                else
                    cells(:,2*p)   = (/ op+p, op+ndx*ndy+j, op+nx_global*ndy+ndx*ndy /)
                    cells(:,2*p+1) = (/ op+p, op+nx_global*ndy+ndx*ndy, op+nx_global*ndy+i*ndy /)
                endif
            enddo
        enddo
    endif


    do i = 0, lc-1
        C(i) = 1.0*DBLE(rank)
    enddo

    schema_version = "1.1"
    dimemsions = "nx_global,ny_global"

    call adios_init_noxml (comm, adios_err)
    call adios_set_max_buffer_size (10) 

    !!
    !! Define the Data group
    !!
    call adios_declare_group (data_group, "tri2d", "", 1, adios_err)
    call adios_select_method (data_group, "MPI", "", "", adios_err)

    ! This example doesn't use varid during writing.
    ! So we simply put 'varid' everywhere.
    ! define a integer
    !call adios_define_var (data_group, "nx_global" &
    !            ,"", adios_integer &
    !            ,"","","", varid)
    !call adios_define_var (data_group, "ny_global" &
    !            ,"", adios_integer &
    !            ,"" ,"","", varid)
    !call adios_define_var (data_group, "nproc" &
    !            ,"", adios_integer &
    !            ,"" ,"","", varid)
    call adios_define_var (data_group, "npoints" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (data_group, "num_cells" &
                ,"", adios_integer &
                ,"","","", varid)
    !call adios_define_var (data_group, "offs_x" &
    !            ,"", adios_integer &
    !            ,"","","", varid)
    !call adios_define_var (data_group, "offs_y" &
    !            ,"", adios_integer &
    !            ,"","" ,"" , varid)
    !call adios_define_var (data_group, "nx_local" &
    !            ,"", adios_integer &
    !            ,"","" , "", varid)
    !call adios_define_var (data_group, "ny_local" &
    !            ,"", adios_integer &
    !            ,"","","", varid)
    call adios_define_var (data_group, "op" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (data_group, "lp" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (data_group, "oc" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (data_group, "lc" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (data_group, "N" &
                    ,"", adios_double &
                    ,"lp", "npoints", "op", varid)
    call adios_define_var (data_group, "C" &
                    ,"", adios_double &
                    ,"lc", "num_cells", "oc", varid)

    !call adios_define_attribute (data_group, "description", "/nproc", adios_string, "Number of writers", "", adios_err)
    call adios_define_attribute (data_group, "description", "/npoints", adios_string, "Number of points", "", adios_err)
    call adios_define_attribute (data_group, "description", "/num_cells", adios_string, "Number of triangles", "", adios_err)

    call adios_define_schema_version (data_group, schema_version)
    call adios_define_mesh_file (data_group, "trimesh", "trimesh.bp")

    call adios_define_var_mesh (data_group, "N", "trimesh")
    call adios_define_var_centering (data_group, "N", "point")
    call adios_define_attribute (data_group, "description", "/N", adios_string, "Node centered data", "", adios_err)
    call adios_define_var_mesh (data_group, "C", "trimesh")
    call adios_define_var_centering (data_group, "C", "cell")
    call adios_define_attribute (data_group, "description", "/C", adios_string, "Cell centered data", "", adios_err)


    !!
    !! Define the Mesh group
    !!
    call adios_declare_group (mesh_group, "trimesh", "", 1, adios_err)
    call adios_select_method (mesh_group, "MPI", "", "", adios_err)
    call adios_define_var (mesh_group, "npoints" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (mesh_group, "num_cells" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (mesh_group, "op" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (mesh_group, "lp" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (mesh_group, "oc" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (mesh_group, "lc" &
                ,"", adios_integer &
                ,"","","", varid)
    call adios_define_var (mesh_group, "points" &
                    ,"", adios_double &
                    ,"2,lp", "2,npoints", "0,op", varid)
    call adios_define_var (mesh_group, "cells" &
                    ,"", adios_integer &
                    ,"3,lc", "3,num_cells", "0,oc", varid)

    call adios_define_attribute (mesh_group, "description", "/npoints", adios_string, "Number of points", "", adios_err)
    call adios_define_attribute (mesh_group, "description", "/num_cells", adios_string, "Number of triangles", "", adios_err)

    call adios_define_schema_version (mesh_group, schema_version)
    call adios_define_mesh_timevarying ("no", mesh_group, "trimesh")
    call adios_define_mesh_unstructured ("points", "cells", "num_cells", "triangle", "", "2", mesh_group, "trimesh")


    !!
    !! Write Data group
    !!
    call adios_open (data_file, "tri2d", "tri2d_f_noxml_sep.bp", "w", comm, adios_err)

    adios_groupsize = 6*4 &   ! int
    + 8 * lp  &              ! double
    + 8 * lc                ! double

    call adios_group_size (data_file, adios_groupsize, adios_totalsize, adios_err)
    
    !call adios_write (data_file, "nproc", nproc, adios_err)
    call adios_write (data_file, "npoints", npoints, adios_err)
    call adios_write (data_file, "num_cells", num_cells, adios_err)
    !call adios_write (data_file, "nx_global", nx_global, adios_err)
    !call adios_write (data_file, "ny_global", ny_global, adios_err)
    !call adios_write (data_file, "offs_x", offs_x, adios_err)
    !call adios_write (data_file, "offs_y", offs_y, adios_err)
    !call adios_write (data_file, "nx_local", nx_local, adios_err)
    !call adios_write (data_file, "ny_local", ny_local, adios_err)
    call adios_write (data_file, "lp", lp, adios_err)
    call adios_write (data_file, "op", op, adios_err)
    call adios_write (data_file, "lc", lc, adios_err)
    call adios_write (data_file, "oc", oc, adios_err)
    call adios_write (data_file, "N", N, adios_err)
    call adios_write (data_file, "C", C, adios_err)

    call adios_close (data_file, adios_err)

    !!
    !! Write Mesh group
    !!
    call adios_open (mesh_file, "trimesh", "trimesh.bp", "w", comm, adios_err)

    adios_groupsize = 6*4 &   ! int
    + 8 * lp * 2 &           ! double
    + 4 * lc * 3 

    call adios_group_size (mesh_file, adios_groupsize, adios_totalsize, adios_err)
    call adios_write (mesh_file, "npoints", npoints, adios_err)
    call adios_write (mesh_file, "num_cells", num_cells, adios_err)
    call adios_write (mesh_file, "lp", lp, adios_err)
    call adios_write (mesh_file, "op", op, adios_err)
    call adios_write (mesh_file, "lc", lc, adios_err)
    call adios_write (mesh_file, "oc", oc, adios_err)
    call adios_write (mesh_file, "points", points, adios_err)
    call adios_write (mesh_file, "cells", cells, adios_err)

    call adios_close (mesh_file, adios_err)


    call MPI_Barrier (comm, ierr)

!    call free (N)
!    call free (points)
!    call free (C)
!    call free (cells)

    call adios_finalize (rank, adios_err)

    call MPI_Finalize (ierr)

end program 

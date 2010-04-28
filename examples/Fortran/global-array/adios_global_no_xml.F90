!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

! ADIOS Fortran Example: write a global array from N processors with gwrite
!
! How to run: mpirun -np <N> adios_global
! Output: adios_global.bp
! ADIOS config file: adios_global.xml
!

program adios_global 
    implicit none
    include 'mpif.h'
    character(len=256)      :: filename = "adios_global.bp"
    integer                 :: rank, size, i, ierr
    integer,parameter       :: NX=10
    integer                 :: O, G
    real*8, dimension(NX)   :: t
    integer                 :: comm

    ! ADIOS variables declarations for matching gwrite_temperature.fh 
    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle
    integer*8               :: m_adios_group

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)

    G = 2 * NX * size
    O = rank * 2 * NX;

    do i = 1, NX
        t(i)  = 10.0*rank+i
    enddo

    call adios_init_noxml (adios_err)
    call adios_allocate_buffer (1, 10, adios_err)

    call adios_declare_group (m_adios_group, "restart", "iter", adios_err)
    call adios_select_method (m_adios_group, "MPI", "", "", adios_err)

    ! define a string
    call adios_define_var (m_adios_group, "Group Name" &
                          ,"", 9 &
                          ,0, 0, 0, adios_err)
    ! define a integer
    call adios_define_var (m_adios_group, "NX" &
                          ,"", 2 &
                          ,0, 0, 0, adios_err)
    ! define a integer
    call adios_define_var (m_adios_group, "G" &
                          ,"", 2 &
                          ,0, 0, 0, adios_err)
    ! define a integer
    call adios_define_var (m_adios_group, "O" &
                          ,"", 2 &
                          ,0, 0, 0, adios_err)
    ! define a global array
    call adios_define_var (m_adios_group, "temperature" &
                          ,"", 6 &
                          ,"NX", "G", "O", adios_err)
    ! define a string attribute
    call adios_define_attribute (m_adios_group, "attr" &
                                ,"", 9, "Oct 1st", 0, adios_err)

    call adios_open (adios_handle, "restart", filename, "w", comm, adios_err)

    adios_groupsize = 14 + 4 + 4 + 4 + NX * 8
    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

    call adios_write (adios_handle, "Group Name", "Restart Group", adios_err)
    call adios_write (adios_handle, "NX", NX, adios_err)
    call adios_write (adios_handle, "G", G, adios_err)
    call adios_write (adios_handle, "O", O, adios_err)
    call adios_write (adios_handle, "temperature", t, adios_err)

    call adios_close (adios_handle, adios_err)

    call MPI_Barrier (comm, ierr)

    call adios_finalize (rank, adios_err)

    call MPI_Finalize ()
end program

!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

! ADIOS Fortran Example: write a global array from N processors using No-XML API. 
! The example also shows how to do multi-writes in ADIOS.
! How to run: mpirun -np <N> global_array_write_noxml_F
! Output: global_array_noxml_F.bp
! ADIOS config file: None
!

program adios_global 
    use adios_write_mod
    implicit none
    include 'mpif.h'
    character(len=256)      :: filename = "global_array_noxml_F.bp"
    integer                 :: rank, size, i, ierr
    integer,parameter       :: NX=10
    integer                 :: O, G
    real*8, dimension(NX)   :: t
    integer                 :: comm

    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle
    integer*8               :: m_adios_group
    integer*8               :: varid

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)

    call adios_init_noxml (comm, adios_err)
    call adios_set_max_buffer_size (10) 

    call adios_declare_group (m_adios_group, "restart", "iter", 1, adios_err)
    call adios_select_method (m_adios_group, "MPI", "", "", adios_err)

    ! This example doesn't use varid during writing.
    ! So we simply put 'varid' everywhere.
    ! define a integer
    call adios_define_var (m_adios_group, "NX" &
                          ,"", 2 &
                          ,"", "", "", varid)
    ! define a integer
    call adios_define_var (m_adios_group, "G" &
                          ,"", 2 &
                          ,"", "", "", varid)
    ! define a integer
    call adios_define_var (m_adios_group, "O" &
                          ,"", 2 &
                          ,"", "", "", varid)
    ! define a global array
    call adios_define_var (m_adios_group, "temperature" &
                          ,"", 6 &
                          ,"NX", "G", "O", varid)

    ! define a integer
    call adios_define_var (m_adios_group, "NX" &
                          ,"", 2 &
                          ,"", "", "", varid) 
    ! define a integer
    call adios_define_var (m_adios_group, "G" &
                          ,"", 2 &
                          ,"", "", "", varid)
    ! define a integer
    call adios_define_var (m_adios_group, "O" &
                          ,"", 2 &
                          ,"", "", "", varid)
    ! define a global array
    call adios_define_var (m_adios_group, "temperature" &
                          ,"", 6 &
                          ,"NX", "G", "O", varid) 

    call adios_open (adios_handle, "restart", filename, "w", comm, adios_err)

    adios_groupsize = 4 + 4 + 4 + NX * 8 &
                    + 4 + 4 + 4 + NX * 8
    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

    G = 2 * NX * size
    O = 2 * NX * rank
    do i = 1, NX
        t(i)  = O + i - 1
    enddo

    call adios_write (adios_handle, "NX", NX, adios_err)
    call adios_write (adios_handle, "G", G, adios_err)
    call adios_write (adios_handle, "O", O, adios_err)
    call adios_write (adios_handle, "temperature", t, adios_err)


    O = 2 * NX * rank + NX
    do i = 1, NX
        t(i)  = O + i - 1
    enddo

    call adios_write (adios_handle, "NX", NX, adios_err)
    call adios_write (adios_handle, "G", G, adios_err)
    call adios_write (adios_handle, "O", O, adios_err)
    call adios_write (adios_handle, "temperature", t, adios_err)

    call adios_close (adios_handle, adios_err)

    call MPI_Barrier (comm, ierr)

    call adios_finalize (rank, adios_err)

    call MPI_Finalize (ierr)
end program

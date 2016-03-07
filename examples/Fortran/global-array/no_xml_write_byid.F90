!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

! ADIOS Fortran Example: write a global array from N processors using No-XML API. 
! The example also shows how to do multi-writes in ADIOS.
! How to run: mpirun -np <N> adios_global_no_xml
! Output: adios_global_no_xml.bp
! ADIOS config file: None
!

program no_xml_write_byid
    use adios_write_mod
    implicit none
    include 'mpif.h'
    character(len=256)      :: filename = "no_xml_write_byid.bp"
    integer                 :: rank, size, i, ierr
    integer,parameter       :: NX=10
    integer                 :: O, G
    real*8, dimension(NX)   :: t
    integer                 :: comm

    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle
    integer*8               :: m_adios_group
    integer*8               :: var_id1, var_id2
    character(len=32)       :: local, global, offset
    !! attributes (from C variables)
    integer, dimension(5)   :: someints = (/ 5,4,3,2,1 /)
    real*8, dimension(5)    :: somedoubles = (/ 5.55555, 4.4444, 3.333, 2.22, 1.1 /)
    character(len=5), dimension(3) :: three_strings = (/ "X    ", "Yy   ", "ZzZ  " /)


    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)

    call adios_init_noxml (comm, adios_err)
    call adios_set_max_buffer_size (10) 

    call adios_declare_group (m_adios_group, "restart", "iter", 1, adios_err)
    call adios_select_method (m_adios_group, "MPI", "", "", adios_err)

    G = 2 * NX * size
    O = 2 * NX * rank

    write (local, "(I2)") NX
    write (global, "(I3)") G
    write (offset, "(I3)") O

    ! define a global array
    call adios_define_var (m_adios_group, "temperature" &
                          ,"", 6 &
                          ,local, global, offset, var_id1)

    call adios_set_transform (var_id1, "identity", adios_err)


    write (offset, "(I3)") O + NX

    ! define a global array
    call adios_define_var (m_adios_group, "temperature" &
                          ,"", 6 &
                          ,local, global, offset, var_id2)

    call adios_set_transform (var_id2, "identity", adios_err)


    !! add some attributes
    call adios_define_attribute_byvalue (m_adios_group, &
            "single_string","", 1, "A single string attribute", adios_err)
    call adios_define_attribute_byvalue (m_adios_group, &
            "three_strings","", 3, three_strings, adios_err)
    call adios_define_attribute_byvalue (m_adios_group, &
            "single_int",   "", 1, someints, adios_err)
    call adios_define_attribute_byvalue (m_adios_group, &
            "single_double","", 1, somedoubles, adios_err)
    call adios_define_attribute_byvalue (m_adios_group, &
            "five_ints",    "", 5, someints, adios_err)
    call adios_define_attribute_byvalue (m_adios_group, &
            "five_double",  "", 5, somedoubles, adios_err)


    call adios_open (adios_handle, "restart", filename, "w", comm, adios_err)

    adios_groupsize =  NX * 8 &
                    +  NX * 8
    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

    do i = 1, NX
        t(i)  = O + i - 1
    enddo

    call adios_write_byid (adios_handle, var_id1, t, adios_err)

    O = 2 * NX * rank + NX
    do i = 1, NX
        t(i)  = O + i - 1
    enddo

    call adios_write_byid (adios_handle, var_id2, t, adios_err)

    call adios_close (adios_handle, adios_err)

    call MPI_Barrier (comm, ierr)

    call adios_finalize (rank, adios_err)

    call MPI_Finalize (ierr)
end program

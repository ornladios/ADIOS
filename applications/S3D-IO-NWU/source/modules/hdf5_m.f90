!
! Author: Wei-keng Liao
!         wkliao@ece.northwestern.edu
!

!==============================================================================
! module for Read-Write Restart files using HDF5
module hdf5_m

#ifdef _BUILD_HDF5_F
    use hdf5
#endif
    implicit none

    PRIVATE
#ifdef _BUILD_HDF5_F
    integer(HSIZE_T) g_sizes(4), subsizes(4), strides(4)
    integer(HSIZE_T) starts(4)
    integer(HID_T)   yspecies_id, yspecies_fspace, yspecies_mspace
    integer(HID_T)   temp_id,     temp_fspace,     temp_mspace
    integer(HID_T)   pressure_id, pressure_fspace, pressure_mspace
    integer(HID_T)   u_id,        u_fspace,        u_mspace
    integer(HID_T)   md_id,       md_fspace
    integer(HID_T)   file_id, plist_id
    character*16     u_name(3)  ! varaible u, velocity variable names
#endif
    public :: hdf5_write, hdf5_read

    contains

!----< handle_err() >-----------------------------------------------------
subroutine handle_err(err_msg, errcode)
    implicit none
    include 'mpif.h'
    integer,       intent(in) :: errcode
    character*(*), intent(in) :: err_msg

    print *, 'Error: ',trim(err_msg)
!   call MPI_Abort(MPI_COMM_SELF, -1);
    return
end subroutine handle_err

#ifdef _BUILD_HDF5_F
!----< h5_common() >------------------------------------------------------
subroutine h5_common
    use topology_m,  only : mypx, mypy, mypz
    use param_m,     only : nx, ny, nz, nx_g, ny_g, nz_g

    implicit none

    ! global array dimensionality
    g_sizes(1) = nx_g
    g_sizes(2) = ny_g
    g_sizes(3) = nz_g

    ! local subarray dimensionality
    subsizes(1) = nx
    subsizes(2) = ny
    subsizes(3) = nz

    ! start offsets of local array in global array
    starts(1) = nx * mypx
    starts(2) = ny * mypy
    starts(3) = nz * mypz
    starts(4) = 0

    ! variable u's names
    u_name(1) = 'u'
    u_name(2) = 'v'
    u_name(3) = 'w'

    return
end subroutine h5_common

!----< hdf5_write_scalar_integer() >--------------------------------------------
subroutine hdf5_write_scalar_integer(name, var)
    use hdf5
    implicit none

    ! declarations passed in
    CHARACTER(LEN=*), INTENT(IN) :: name   ! Name of the dataset 
    INTEGER,          INTENT(IN) :: var

    ! local variables
    integer ierr
    call h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, md_fspace, md_id, ierr)
    ! since buf(time) is a scalar, dims(subsizes) agument is ignored
    call h5dwrite_f(md_id, H5T_NATIVE_INTEGER, var, subsizes, ierr)
    call h5dclose_f(md_id, ierr)

    return
end subroutine hdf5_write_scalar_integer

!----< hdf5_write_scalar_double() >--------------------------------------------
subroutine hdf5_write_scalar_double(name, var)
    use hdf5
    implicit none

    ! declarations passed in
    CHARACTER(LEN=*), INTENT(IN) :: name   ! Name of the dataset 
    DOUBLE PRECISION, INTENT(IN) :: var

    ! local variables
    integer ierr
    call h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, md_fspace, md_id, ierr)
    ! since buf(time) is a scalar, dims(subsizes) agument is ignored
    call h5dwrite_f(md_id, H5T_NATIVE_DOUBLE, var, subsizes, ierr)
    call h5dclose_f(md_id, ierr)

    return
end subroutine hdf5_write_scalar_double

!----< hdf5_write_scalar_string() >--------------------------------------------
subroutine hdf5_write_scalar_string(name, type_id, var)
    use hdf5
    implicit none

    ! declarations passed in
    CHARACTER(LEN=*), INTENT(IN) :: name   ! Name of the dataset 
    CHARACTER(LEN=*), INTENT(IN) :: var
    INTEGER(HID_T),   INTENT(IN) :: type_id  ! Dataspace identifier

    ! local variables
    integer ierr

    call h5dcreate_f(file_id, name, type_id, md_fspace, md_id, ierr)
    call h5dwrite_f(md_id, type_id, var, subsizes, ierr)
    call h5dclose_f(md_id, ierr)

    return
end subroutine hdf5_write_scalar_string
#endif

!----< hdf5_write() >---------------------------------------------------------
subroutine hdf5_write(filename)
    use topology_m,  only : gcomm
    use param_m,     only : nsc, n_elem, n_spec, n_reac, ntr
    use variables_m, only : temp, pressure, yspecies, u
    use runtime_m,   only : time, tstep, time_save
    use bc_m,        only : pout
    use chemkin_m,   only : species_name, n_species, molwt, element_name, n_elements
    use reference_m

#ifdef _BUILD_HDF5_F
    use hdf5
#endif

    implicit none
    include 'mpif.h'

    ! declarations passed in
    character*(*), intent(in) :: filename

#ifdef _BUILD_HDF5_F
    ! local variables
    integer i, ierr, mpi_io_fh, file_info
    INTEGER(HSIZE_T), DIMENSION(1) :: dims
    INTEGER(HID_T) :: type_id      ! Attribute Dataspace identifier
    character*16 str
    character*16, allocatable :: species_name2(:)

    call h5_common

    ! Initialize FORTRAN predefined datatypes --------------------------------
    call h5open_f(ierr)

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call MPI_Info_create(file_info, ierr)
    call MPI_Info_set(file_info, 'romio_ds_write', 'disable', ierr)
!   call MPI_Info_set(file_info, 'romio_no_indep_rw', 'true', ierr)
    call h5pset_fapl_mpio_f(plist_id, gcomm, file_info, ierr)

    ! Create the file collectively -------------------------------------------
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, ierr, &
                     access_prp = plist_id)
    call MPI_Info_free(file_info, ierr)

    ! if (myid .EQ. 0) then
    !     the HDF5 function below is not implemented yet
    !     call h5fget_vfd_handle_f(file_id, plist_id, mpi_io_fh, ierr)
    !     call MPI_File_get_info(mpi_io_fh, mpi_info_used, ierr)
    !     call show_file_info(mpi_info_used)
    !     call MPI_Info_free(mpi_info_used, ierr)
    ! endif

    call h5pclose_f(plist_id, ierr)

    ! Save metadata as HDF5 dataset -----------------------------------
    call h5screate_f(H5S_SCALAR_F, md_fspace, ierr)

    call hdf5_write_scalar_integer('number_of_elements_in_reaction_mechansim', n_elem)
    call hdf5_write_scalar_integer('number_of_species_in_reaction_mechansim', n_spec)
    call hdf5_write_scalar_integer('number_of_steps_in_reaction_mechansims', n_reac/2)
    call hdf5_write_scalar_integer('number_of_reaction_third-body_reactions', ntr)

    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierr)
    call h5tset_size_f(type_id, 16, ierr)

    do i=1,n_elements,1
 1001 format('element_name_',I1)
        write(str, 1001) i
        call hdf5_write_scalar_string(str, type_id, element_name(i))
    enddo

    do i=1,n_species,1
       call hdf5_write_scalar_double('molecular_weight:'//species_name(i), molwt(i))
    enddo

    call hdf5_write_scalar_double('universal_gas_constant',            univ_gascon)
    call hdf5_write_scalar_double('freestream_temperature',            t_o)
    call hdf5_write_scalar_double('reference_ratio_of_specifice_heats',g_ref)
    call hdf5_write_scalar_double('reference_speed_of_sound',          a_ref)
    call hdf5_write_scalar_double('reference_density',                 rho_ref)
    call hdf5_write_scalar_double('reference_conductivity',            lambda_ref)
    call hdf5_write_scalar_double('reference_temperature',             t_ref)
    call hdf5_write_scalar_double('reference_pressure',                p_ref)
    call hdf5_write_scalar_double('standard_atmospheric_pressure',     pres_atm)
    call hdf5_write_scalar_double('reference_time',                    time_ref)
    call hdf5_write_scalar_double('reference_specific_heat',           cp_ref)
    call hdf5_write_scalar_double('reference_length',                  l_ref)
    call hdf5_write_scalar_double('reference_viscosity',               mu_ref)
    call hdf5_write_scalar_double('acoustic_Reynolds_number',          re)
    call hdf5_write_scalar_double('Mach_number',                       mach_no)
    call hdf5_write_scalar_double('convective_Reynolds_number',        re*mach_no)

    call hdf5_write_scalar_double('time', time)
    call hdf5_write_scalar_double('tstep', tstep)
    call hdf5_write_scalar_double('time_save', time_save)
    call hdf5_write_scalar_double('pout', pout)

    call h5sclose_f(md_fspace, ierr)

    ! Create the file data space for all arrays ------------------------------
    g_sizes(4) = nsc+1
    call h5screate_simple_f(4, g_sizes, yspecies_fspace, ierr)
    call h5screate_simple_f(3, g_sizes, temp_fspace,     ierr)
    call h5screate_simple_f(3, g_sizes, pressure_fspace, ierr)
    g_sizes(4) = 3
    call h5screate_simple_f(4, g_sizes, u_fspace,        ierr)

    ! Create the array datasets in the file with default properties -----------
    call h5dcreate_f(file_id, 'yspecies', H5T_NATIVE_DOUBLE, yspecies_fspace, yspecies_id, ierr)
    call h5dcreate_f(file_id, 'temp',     H5T_NATIVE_DOUBLE, temp_fspace,     temp_id,     ierr)
    call h5dcreate_f(file_id, 'pressure', H5T_NATIVE_DOUBLE, pressure_fspace, pressure_id, ierr)
    call h5dcreate_f(file_id, 'u',        H5T_NATIVE_DOUBLE, u_fspace,        u_id,        ierr)

    allocate(species_name2(nsc+1))
    do i=1,nsc+1
        species_name2(i) = 'Y-'//species_name(i)
    enddo

    dims(1)=nsc+1
    call h5screate_simple_f(1, dims, md_fspace, ierr)
    CALL h5acreate_f(yspecies_id, 'specie names', type_id, md_fspace, md_id, ierr)
    call h5awrite_f(md_id, type_id, species_name2, dims, ierr)
    call h5aclose_f(md_id, ierr)
    call h5sclose_f(md_fspace, ierr)
    deallocate(species_name2)

    dims(1)=3
    call h5screate_simple_f(1, dims, md_fspace, ierr)
    CALL h5acreate_f(u_id, 'velocity component', type_id, md_fspace, md_id, ierr)
    call h5awrite_f(md_id, type_id, u_name, dims, ierr)
    call h5aclose_f(md_id, ierr)
    call h5sclose_f(md_fspace, ierr)

    ! Create the memory data space for local subarrays -----------------------
    subsizes(4) = nsc+1
    call h5screate_simple_f(4, subsizes, yspecies_mspace, ierr)
    call h5screate_simple_f(3, subsizes, temp_mspace,     ierr)
    call h5screate_simple_f(3, subsizes, pressure_mspace, ierr)
    subsizes(4) = 3
    call h5screate_simple_f(4, subsizes, u_mspace,        ierr)

    ! Select hyperslab in the file -------------------------------------------
    strides(:) = 1

    subsizes(4) = nsc+1
    call h5sselect_hyperslab_f(yspecies_fspace, H5S_SELECT_SET_F, starts, &
                               subsizes, ierr, stride = strides)
    call h5sselect_hyperslab_f(temp_fspace, H5S_SELECT_SET_F, starts, &
                               subsizes, ierr, stride = strides)
    call h5sselect_hyperslab_f(pressure_fspace, H5S_SELECT_SET_F, starts, &
                               subsizes, ierr, stride = strides)
    subsizes(4) = 3
    call h5sselect_hyperslab_f(u_fspace, H5S_SELECT_SET_F, starts, &
                               subsizes, ierr, stride = strides)

    ! Create property list for collective dataset write ----------------------
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Write the dataset yspecies collectively.
    subsizes(4) = nsc+1
    call h5dwrite_f(yspecies_id, H5T_NATIVE_DOUBLE, yspecies, subsizes, &
                    ierr, yspecies_mspace, yspecies_fspace, plist_id)

    ! Write the dataset temp collectively.
    call h5dwrite_f(temp_id, H5T_NATIVE_DOUBLE, temp, subsizes, ierr, &
                    temp_mspace, temp_fspace, plist_id)

    ! Write the dataset pressure collectively.
    call h5dwrite_f(pressure_id, H5T_NATIVE_DOUBLE, pressure, subsizes, &
                    ierr, pressure_mspace, pressure_fspace, plist_id)

    ! Write the dataset u collectively.
    subsizes(4) = 3
    call h5dwrite_f(u_id, H5T_NATIVE_DOUBLE, u, subsizes, ierr, &
                    u_mspace, u_fspace, plist_id)

    ! Close dataspaces
    call h5sclose_f(yspecies_fspace, ierr)
    call h5sclose_f(yspecies_mspace, ierr)
    call h5sclose_f(temp_fspace,     ierr)
    call h5sclose_f(temp_mspace,     ierr)
    call h5sclose_f(pressure_fspace, ierr)
    call h5sclose_f(pressure_mspace, ierr)
    call h5sclose_f(u_fspace,        ierr)
    call h5sclose_f(u_mspace,        ierr)

    ! Close the datasets
    call h5dclose_f(yspecies_id, ierr)
    call h5dclose_f(temp_id,     ierr)
    call h5dclose_f(pressure_id, ierr)
    call h5dclose_f(u_id,        ierr)

    ! Close the property list.
    call h5pclose_f(plist_id, ierr)

    ! Close the file.
    call h5fclose_f(file_id, ierr)

    ! Close FORTRAN predefined datatypes.
    call h5close_f(ierr)
#endif

    return
end subroutine hdf5_write

!----< hdf5_read() >---------------------------------------------------------
subroutine hdf5_read(filename)
    use topology_m,  only : gcomm
    use param_m,     only : nsc
    use variables_m, only : temp, pressure, yspecies, u
    use runtime_m,   only : time, tstep, time_save
    use bc_m,        only : pout

#ifdef _BUILD_HDF5_F
    use hdf5
#endif

    implicit none
    include 'mpif.h'

    ! declarations passed in
    character*(*), intent(in) :: filename

#ifdef _BUILD_HDF5_F
    ! local variables
    integer i, ierr, file_info

    call h5_common

    ! Initialize FORTRAN predefined datatypes --------------------------------
    call h5open_f(ierr)

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call MPI_Info_create(file_info, ierr)
    call MPI_Info_set(file_info, 'romio_ds_write', 'disable', ierr)
!   call MPI_Info_set(file_info, 'romio_no_indep_rw', 'true', ierr)
    call h5pset_fapl_mpio_f(plist_id, gcomm, file_info, ierr)
    call MPI_Info_free(file_info, ierr)

    ! Open the file collectively ---------------------------------------------
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, ierr, &
                   access_prp = plist_id)

    call h5pclose_f(plist_id, ierr)

    ! Get timing metadata as HDF5 dataset -----------------------------------
    ! since all metadata attributes are scalar, dims(subsizes) argument is ignored
    call h5dopen_f(file_id, 'time', md_id, ierr)
    call h5dread_f(md_id, H5T_NATIVE_DOUBLE, time, subsizes, ierr)
    call h5dclose_f(md_id, ierr)

    call h5dopen_f(file_id, 'tstep', md_id, ierr)
    call h5dread_f(md_id, H5T_NATIVE_DOUBLE, tstep, subsizes, ierr)
    call h5dclose_f(md_id, ierr)

    call h5dopen_f(file_id, 'time_save', md_id, ierr)
    call h5dread_f(md_id, H5T_NATIVE_DOUBLE, time_save, subsizes, ierr)
    call h5dclose_f(md_id, ierr)

    call h5dopen_f(file_id, 'pout', md_id, ierr)
    call h5dread_f(md_id, H5T_NATIVE_DOUBLE, pout, subsizes, ierr)
    call h5dclose_f(md_id, ierr)

    ! open all datasets ------------------------------------------------------
    call h5dopen_f(file_id, 'yspecies', yspecies_id, ierr)
    call h5dopen_f(file_id, 'temp',     temp_id,     ierr)
    call h5dopen_f(file_id, 'pressure', pressure_id, ierr)
    call h5dopen_f(file_id, 'u',        u_id,        ierr)

    ! get file spaces for all datasets ---------------------------------------
    call h5dget_space_f(yspecies_id, yspecies_fspace, ierr)
    call h5dget_space_f(temp_id,     temp_fspace,     ierr)
    call h5dget_space_f(pressure_id, pressure_fspace, ierr)
    call h5dget_space_f(u_id,        u_fspace,        ierr)

    ! Create the memory data space for local subarrays -----------------------
    subsizes(4) = nsc+1
    call h5screate_simple_f(4, subsizes, yspecies_mspace, ierr)
    call h5screate_simple_f(3, subsizes, temp_mspace,     ierr)
    call h5screate_simple_f(3, subsizes, pressure_mspace, ierr)
    subsizes(4) = 3
    call h5screate_simple_f(4, subsizes, u_mspace,        ierr)

    ! Select hyperslab in the file ------------------------------------------
    strides(:) = 1

    subsizes(4) = nsc+1
    call h5sselect_hyperslab_f(yspecies_fspace, H5S_SELECT_SET_F, starts, &
                               subsizes, ierr, stride = strides)

    call h5sselect_hyperslab_f(temp_fspace, H5S_SELECT_SET_F, starts, &
                               subsizes, ierr, stride = strides)

    call h5sselect_hyperslab_f(pressure_fspace, H5S_SELECT_SET_F, starts, &
                               subsizes, ierr, stride = strides)
    subsizes(4) = 3
    call h5sselect_hyperslab_f(u_fspace, H5S_SELECT_SET_F, starts, &
                               subsizes, ierr, stride = strides)

    ! Create property list for collective dataset read -----------------------
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Read the dataset yspecies collectively.
    subsizes(4) = nsc+1
    call h5dread_f(yspecies_id, H5T_NATIVE_DOUBLE, yspecies, subsizes, &
                   ierr, yspecies_mspace, yspecies_fspace, plist_id)

    ! Read the dataset temp collectively.
    call h5dread_f(temp_id, H5T_NATIVE_DOUBLE, temp, subsizes, ierr, &
                   temp_mspace, temp_fspace, plist_id)

    ! Read the dataset pressure collectively.
    call h5dread_f(pressure_id, H5T_NATIVE_DOUBLE, pressure, subsizes, &
                   ierr, pressure_mspace, pressure_fspace, plist_id)

    ! Read the dataset u collectively.
    subsizes(4) = 3
    call h5dread_f(u_id, H5T_NATIVE_DOUBLE, u, subsizes, ierr, &
                   u_mspace, u_fspace, plist_id)

    ! Close dataspaces
    call h5sclose_f(yspecies_fspace, ierr)
    call h5sclose_f(yspecies_mspace, ierr)
    call h5sclose_f(temp_fspace,     ierr)
    call h5sclose_f(temp_mspace,     ierr)
    call h5sclose_f(pressure_fspace, ierr)
    call h5sclose_f(pressure_mspace, ierr)
    call h5sclose_f(u_fspace, ierr)
    call h5sclose_f(u_mspace, ierr)

    ! Close the datasets
    call h5dclose_f(yspecies_id, ierr)
    call h5dclose_f(temp_id,     ierr)
    call h5dclose_f(pressure_id, ierr)
    call h5dclose_f(u_id,        ierr)

    ! Close the property list.
    call h5pclose_f(plist_id, ierr)

    ! Close the file.
    call h5fclose_f(file_id, ierr)

    ! Close FORTRAN predefined datatypes.
    call h5close_f(ierr)

#endif

    return
end subroutine hdf5_read

!==============================================================================
end module hdf5_m


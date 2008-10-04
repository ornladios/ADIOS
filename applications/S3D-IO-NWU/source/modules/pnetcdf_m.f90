!
! Author: Wei-keng Liao
!         wkliao@ece.northwestern.edu
!

!==============================================================================
module pnetcdf_m
    ! module for Read-Write Restart files using Parallel NetCDF

    implicit none
    include 'mpif.h'
#ifdef _BUILD_PNETCDF_F
#   include "pnetcdf.inc"
#endif

    PRIVATE
    integer(MPI_OFFSET_KIND) g_sizes(4), subsizes(4), starts(4)
    integer dimids(4)
    integer yspecies_id, pressure_id, temp_id, u_id

    public :: pnetcdf_write, pnetcdf_read

    contains

!----< handle_err() >---------------------------------------------------------
subroutine handle_err(err_msg, errcode)
    implicit none

    integer,       intent(in) :: errcode
    character*(*), intent(in) :: err_msg
     
#ifdef _BUILD_PNETCDF_F
    print *, 'Error: ',trim(err_msg),' ',nfmpi_strerror(errcode)
#endif
    call MPI_Abort(MPI_COMM_SELF, -1);
    return
end subroutine handle_err

!----< pnetcdf_common() >------------------------------------------------------
subroutine pnetcdf_common
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
    ! note that Fortran array index starts with 1
    starts(1) = nx * mypx + 1
    starts(2) = ny * mypy + 1
    starts(3) = nz * mypz + 1
    starts(4) = 1

    return
end subroutine pnetcdf_common

!----< pnetcdf_write() >------------------------------------------------------
subroutine pnetcdf_write(filename)
    use topology_m,  only : gcomm
    use param_m,     only : nsc, n_elem, n_spec, n_reac, ntr
    use variables_m, only : temp, pressure, yspecies, u
    use runtime_m,   only : time, tstep, time_save
    use bc_m,        only : pout
    use chemkin_m,   only : species_name, n_species, molwt, element_name, n_elements
    use reference_m

    implicit none

    ! declarations passed in
    character*(*), intent(in) :: filename

#ifdef _BUILD_PNETCDF_F
    ! local variables
    integer i, ierr, cmode, file_info, ncid
    character*32 str

    call pnetcdf_common

    ! set up some MPI I/O hints for further enhancement
    call MPI_Info_create(file_info, ierr)
    call MPI_Info_set(file_info, 'romio_ds_write', 'disable', ierr)
!   call MPI_Info_set(file_info, 'romio_no_indep_rw', 'true', ierr)

    ! create file and pass in the MPI hint
    cmode = NF_CLOBBER + NF_64BIT_OFFSET
    ierr = nfmpi_create(gcomm, trim(filename), cmode, file_info, ncid)
    if (ierr .ne. 0) call handle_err('nfmpi_create', ierr)
    call MPI_Info_free(file_info, ierr)

    ! Save metadata as global attributes ---------------------------------
    ierr = nfmpi_put_att_int(ncid, NF_GLOBAL, 'number_of_elements_in_reaction_mechansim', NF_INT, 1, n_elem)
    ierr = nfmpi_put_att_int(ncid, NF_GLOBAL, 'number_of_species_in_reaction_mechansim', NF_INT, 1, n_spec)
    ierr = nfmpi_put_att_int(ncid, NF_GLOBAL, 'number_of_steps_in_reaction_mechansims', NF_INT, 1, n_reac/2)
    ierr = nfmpi_put_att_int(ncid, NF_GLOBAL, 'number_of_reaction_third-body_reactions', NF_INT, 1, ntr)

    do i=1, n_elements
 1001 format('element_name_',I1)
        write(str, 1001) i
        ierr = nfmpi_put_att_text(ncid, NF_GLOBAL, str, 16, element_name(i))
    enddo

    do i=1,n_species,1
       ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'molecular_weight:'//species_name(i), NF_DOUBLE, 1, molwt(i))
    enddo

    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'universal_gas_constant',            NF_DOUBLE, 1, univ_gascon)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'freestream_temperature',            NF_DOUBLE, 1, t_o)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_ratio_of_specifice_heats',NF_DOUBLE, 1, g_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_speed_of_sound',          NF_DOUBLE, 1, a_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_density',                 NF_DOUBLE, 1, rho_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_conductivity',            NF_DOUBLE, 1, lambda_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_temperature',             NF_DOUBLE, 1, t_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_pressure',                NF_DOUBLE, 1, p_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'standard_atmospheric_pressure',     NF_DOUBLE, 1, pres_atm)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_time',                    NF_DOUBLE, 1, time_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_specific_heat',           NF_DOUBLE, 1, cp_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_length',                  NF_DOUBLE, 1, l_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'reference_viscosity',               NF_DOUBLE, 1, mu_ref)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'acoustic_Reynolds_number',          NF_DOUBLE, 1, re)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'Mach_number',                       NF_DOUBLE, 1, mach_no)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'convective_Reynolds_number',        NF_DOUBLE, 1, re*mach_no)

    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'time',      NF_DOUBLE, 1, time)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'tstep',     NF_DOUBLE, 1, tstep)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'time_save', NF_DOUBLE, 1, time_save)
    ierr = nfmpi_put_att_double(ncid, NF_GLOBAL, 'pout',      NF_DOUBLE, 1, pout)

    ! define X-Y-Z dimensions of the global array
    ierr = nfmpi_def_dim(ncid, 'nx_g',   g_sizes(1), dimids(1))
    if (ierr .ne. 0) call handle_err('nfmpi_def_dim on x', ierr)
    ierr = nfmpi_def_dim(ncid, 'ny_g',   g_sizes(2), dimids(2))
    if (ierr .ne. 0) call handle_err('nfmpi_def_dim on y', ierr)
    ierr = nfmpi_def_dim(ncid, 'nz_g',   g_sizes(3), dimids(3))
    if (ierr .ne. 0) call handle_err('nfmpi_def_dim on z', ierr)

    ! define variable yspecies
    g_sizes(4) = nsc + 1;
    ierr = nfmpi_def_dim(ncid, 'number_of_species', g_sizes(4), dimids(4))
    ierr = nfmpi_def_var(ncid, 'yspecies', NF_DOUBLE, 4, dimids, yspecies_id)
    do i=1, nsc+1
 1002 format('specie_name_',I2.2)
        write(str, 1002) i
        ierr = nfmpi_put_att_text(ncid, yspecies_id, str, 18, 'Y-'//species_name(i))
    enddo

    ! define variable temp
    ierr = nfmpi_def_var(ncid, 'temp', NF_DOUBLE, 3, dimids, temp_id)
    if (ierr .ne. 0) call handle_err('nfmpi_def_var on temp', ierr)

    ! define variable pressure
    ierr = nfmpi_def_var(ncid, 'pressure', NF_DOUBLE, 3, dimids, pressure_id)
    if (ierr .ne. 0) call handle_err('nfmpi_def_var on pressure', ierr)

    ! define variable u
    g_sizes(4) = 3
    ierr = nfmpi_def_dim(ncid, 'number_of_velocity_components', g_sizes(4), dimids(4))
    ierr = nfmpi_def_var(ncid, 'u', NF_DOUBLE, 4, dimids, u_id)
    ierr = nfmpi_put_att_text(ncid, u_id, 'velocity_component_1', 1, 'u')
    ierr = nfmpi_put_att_text(ncid, u_id, 'velocity_component_2', 1, 'v')
    ierr = nfmpi_put_att_text(ncid, u_id, 'velocity_component_3', 1, 'w')

    ! end of define mode
    ierr = nfmpi_enddef(ncid)
    if (ierr .ne. 0) call handle_err('nfmpi_enddef', ierr)

    !---- write array yspecies
    subsizes(4) = nsc + 1
    ierr = nfmpi_put_vara_double_all(ncid, yspecies_id, starts, subsizes, yspecies)
    if (ierr .ne. 0) call handle_err('nfmpi_put_vara_double_all on yspecies', ierr)

    !---- write array temp
    ierr = nfmpi_put_vara_double_all(ncid, temp_id, starts, subsizes, temp)
    if (ierr .ne. 0) call handle_err('nfmpi_put_vara_double_all on temp', ierr)

    !---- write array pressure
    ierr = nfmpi_put_vara_double_all(ncid, pressure_id, starts, subsizes, pressure)
    if (ierr .ne. 0) call handle_err('nfmpi_put_vara_double_all on pressure', ierr)

    !---- write array u
    subsizes(4) = 3
    ierr = nfmpi_put_vara_double_all(ncid, u_id, starts, subsizes, u)
    if (ierr .ne. 0) call handle_err('nfmpi_put_vara_double_all on u', ierr)

    ierr = nfmpi_close(ncid)
    if (ierr .ne. 0) call handle_err('nfmpi_close', ierr)

#endif

    return
end subroutine pnetcdf_write

!----< pnetcdf_read() >---------------------------------------------------------
subroutine pnetcdf_read(filename)
    use topology_m,  only : gcomm
    use param_m,     only : nsc
    use variables_m, only : temp, pressure, yspecies, u
    use runtime_m,   only : time, tstep, time_save
    use bc_m,        only : pout
    use chemkin_m,   only : species_name

    implicit none

    ! declarations passed in
    character*(*), intent(in) :: filename

#ifdef _BUILD_PNETCDF_F
    ! local variables
    integer i, ierr, cmode, file_info, ncid

    call pnetcdf_common

    ! set up some MPI I/O hints for further enhancement
    call MPI_Info_create(file_info, ierr)
    call MPI_Info_set(file_info, 'romio_ds_write', 'disable', ierr)
!   call MPI_Info_set(file_info, 'romio_no_indep_rw', 'true', ierr)

    ! open file and pass in the MPI hint
    cmode = NF_NOWRITE + NF_64BIT_OFFSET
    ierr = nfmpi_open(gcomm, trim(filename), cmode, file_info, ncid)
    if (ierr .ne. 0) call handle_err('nfmpi_open', ierr)
    call MPI_Info_free(file_info, ierr)

    ! Get timing metadata global attributes ---------------------------------
    ierr = nfmpi_get_att_double(ncid, NF_GLOBAL, 'time', time)
    ierr = nfmpi_get_att_double(ncid, NF_GLOBAL, 'tstep', tstep)
    ierr = nfmpi_get_att_double(ncid, NF_GLOBAL, 'time_save', time_save)
    ierr = nfmpi_get_att_double(ncid, NF_GLOBAL, 'pout', pout)

    ! inquire variable id for yspecies
    ierr = nfmpi_inq_varid(ncid, 'yspecies', yspecies_id)
    if (ierr .ne. 0) call handle_err('nfmpi_inq_varid on yspecies', ierr)

    ! inquire variable temp id
    ierr = nfmpi_inq_varid(ncid, 'temp', temp_id)
    if (ierr .ne. 0) call handle_err('nfmpi_inq_varid on temp', ierr)

    ! inquire variable pressure id
    ierr = nfmpi_inq_varid(ncid, 'pressure', pressure_id)
    if (ierr .ne. 0) call handle_err('nfmpi_inq_varid on pressure', ierr)

    ! inquire variable u id
    ierr = nfmpi_inq_varid(ncid, 'u', u_id)
    if (ierr .ne. 0) call handle_err('nfmpi_inq_varid on u', ierr)

    !---- read array yspecies
    subsizes(4) = nsc+1
    ierr = nfmpi_get_vara_double_all(ncid, yspecies_id, starts, subsizes, yspecies)
    if (ierr .ne. 0) call handle_err('nfmpi_get_vara_double_all on yspecies', ierr)

    !---- read array temp
    ierr = nfmpi_get_vara_double_all(ncid, temp_id, starts, subsizes, temp)
    if (ierr .ne. 0) call handle_err('nfmpi_get_vara_double_all on temp', ierr)

    !---- read array pressure
    ierr = nfmpi_get_vara_double_all(ncid, pressure_id, starts, subsizes, pressure)
    if (ierr .ne. 0) call handle_err('nfmpi_get_vara_double_all on pressure', ierr)

    !---- read array u
    subsizes(4) = 3
    ierr = nfmpi_get_vara_double_all(ncid, u_id, starts, subsizes, u)
    if (ierr .ne. 0) call handle_err('nfmpi_get_vara_double_all on u', ierr)

    ierr = nfmpi_close(ncid)
    if (ierr .ne. 0) call handle_err('nfmpi_close', ierr)

#endif

    return
end subroutine pnetcdf_read

!==============================================================================
end module pnetcdf_m

!
! Author: Wei-keng Liao
!         wkliao@ece.northwestern.edu
!

!==============================================================================
module mpi_io_m
! module for Read-Write Restart files using MPI I/O

    implicit none

    integer io_method, file_info, threeD_ftype

    private :: file_info, threeD_ftype, handle_err

  contains

!----< handle_err() >----------------------------------------------------------
subroutine handle_err(err_msg, errcode)
    implicit none
    include 'mpif.h'
    integer,       intent(in) :: errcode
    character*(*), intent(in) :: err_msg
    integer errorclass, ierr, resultlen
    character*256  err_str

    call MPI_Error_string(errcode, err_str, resultlen, ierr)
    print *, 'Error: ',trim(err_msg),' error string=',err_str
    call MPI_Abort(MPI_COMM_SELF, -1);
    return
end subroutine handle_err

!----< mpi_io_set_filetype() >------------------------------------------------
subroutine mpi_io_set_filetype(flag)
    use topology_m,  only : npes, mypx, mypy, mypz
    use param_m,     only : nx, ny, nz, nx_g, ny_g, nz_g
    implicit none
    include 'mpif.h'

    ! declarations passed in
    integer, intent(in) :: flag ! 0 for creating data type, -1 for freeing it up

    ! local variables
    integer       g_sizes(3), subsizes(3), starts(3)
    integer       ierr, errorclass
    character*16  int_str
  
    ! free up info and file type
    if (flag .eq. -1) then
        if (file_info .ne. MPI_INFO_NULL) then
            call MPI_Info_free(file_info, ierr)
            file_info = MPI_INFO_NULL
        endif
        call MPI_Type_free(threeD_ftype, ierr)
        return
    endif

    ! global array dimensions
    g_sizes(1) = nx_g
    g_sizes(2) = ny_g
    g_sizes(3) = nz_g

    ! local array dimensions
    subsizes(1) = nx
    subsizes(2) = ny
    subsizes(3) = nz

    ! local array's start offsets in global array
    starts(1) = nx * mypx
    starts(2) = ny * mypy
    starts(3) = nz * mypz

    ! define local 3D sub-array data type
    call MPI_Type_create_subarray(3, g_sizes, subsizes, &
                  starts, MPI_ORDER_FORTRAN, MPI_REAL8, &
                  threeD_ftype, ierr)
    call MPI_Type_commit(threeD_ftype, ierr)
    ! threeD_ftype will be used to define MPI file view

    ! set up some MPI I/O hints for further enhancement
    call MPI_Info_create(file_info, ierr)
    call MPI_Info_set(file_info, 'romio_ds_write', 'disable', ierr)
    call MPI_Info_set(file_info, 'romio_no_indep_rw', 'true', ierr)

    return
end subroutine mpi_io_set_filetype

!----< mpi_file_io() >---------------------------------------------------------
subroutine mpi_file_io(filename,rw)
    use topology_m,  only : gcomm
    use param_m,     only : nx, ny, nz, nsc
    use variables_m, only : temp, pressure, yspecies, u
    implicit none
    include 'mpif.h'

    ! declarations passed in
    character*1,   intent(in)    :: rw    ! 'r' for read, 'w' for write
    character*100, intent(inout) :: filename

    ! local variables
    integer nx_ny_nz
    integer fp, ierr, errorclass, open_mode
    integer mstatus(MPI_STATUS_SIZE)
    integer(MPI_OFFSET_KIND) iOffset

    nx_ny_nz = nx * ny * nz

    open_mode = MPI_MODE_WRONLY+MPI_MODE_CREATE
    if (rw .EQ. 'r') open_mode = MPI_MODE_RDONLY

    ! open a shared file with an MPI hint
    call MPI_File_open(gcomm, trim(filename), open_mode, file_info, fp, ierr)
    if (ierr .ne. MPI_SUCCESS) call handle_err('MPI_File_open',ierr)

    iOffset = 0
    call MPI_File_set_view(fp, iOffset, MPI_REAL8, threeD_ftype, 'native', MPI_INFO_NULL, ierr)
    if (ierr .ne. MPI_SUCCESS) call handle_err('MPI_File_set_view',ierr)

    !---- write/read array yspecies
    if (rw .EQ. 'w') then
        call MPI_File_write_all(fp, yspecies, nx_ny_nz*(nsc+1), MPI_REAL8, mstatus, ierr)
    else
        call MPI_File_read_all(fp, yspecies, nx_ny_nz*(nsc+1), MPI_REAL8, mstatus, ierr)
    endif
    if (ierr .ne. MPI_SUCCESS) call handle_err('MPI write/read array: yspecies',ierr)

    !---- write/read array temp
    if (rw .EQ. 'w') then
        call MPI_File_write_all(fp, temp, nx_ny_nz, MPI_REAL8, mstatus, ierr)
    else
        call MPI_File_read_all(fp, temp, nx_ny_nz, MPI_REAL8, mstatus, ierr)
    endif
    if (ierr .ne. MPI_SUCCESS) call handle_err('MPI write/read array: temp',ierr)

    !---- write/read array pressure
    if (rw .EQ. 'w') then
        call MPI_File_write_all(fp, pressure, nx_ny_nz, MPI_REAL8, mstatus, ierr)
    else
        call MPI_File_read_all(fp, pressure, nx_ny_nz, MPI_REAL8, mstatus, ierr)
    endif
    if (ierr .ne. MPI_SUCCESS) call handle_err('MPI write/read array: pressure',ierr)

    !---- write/read array u
    if (rw .EQ. 'w') then
        call MPI_File_write_all(fp, u, nx_ny_nz*3, MPI_REAL8, mstatus, ierr)
    else
        call MPI_File_read_all(fp, u, nx_ny_nz*3, MPI_REAL8, mstatus, ierr)
    endif
    if (ierr .ne. MPI_SUCCESS) call handle_err('MPI write/read array: u',ierr)

    !---- close file
    call MPI_File_close(fp, ierr)
    if (ierr .ne. MPI_SUCCESS) call handle_err('MPI_File_close',ierr)

    return
end subroutine mpi_file_io

!==============================================================================
end module mpi_io_m

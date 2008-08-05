program adios_test
    implicit none
    include 'mpif.h'
    character (len=200) :: group
    character (len=200) :: filename
    integer :: group_comm
    integer :: ierr
    integer :: rank

    ! write vars
    integer :: small_int
    integer*8 :: big_int
    real :: small_real
    real*8 :: big_real
    integer :: z_size
    real :: z_array (2)

    ! read vars
    integer :: r_small_int
    integer*8 :: r_big_int
    real :: r_small_real
    real*8 :: r_big_real
    integer :: r_z_size
    real :: r_z_array (2)

    group = "restart"
    filename = "restart.bp"
    small_int = 10
    big_int = 4294967296
    small_real = 0.3
    big_real = 0.00000000000004
    z_size = 2
    z_array (1) = 11.1
    z_array (2) = 22.2
    r_z_size = 2

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)

    call adios_init ("config_fortran.xml"//char(0), ierr)

    call test_write (group, filename, group_comm, small_int, big_int, small_real, big_real, z_size, z_array)

    call MPI_Barrier (MPI_COMM_WORLD, ierr)

    call test_read (group, filename, group_comm, r_small_int, r_big_int, r_small_real, r_big_real, r_z_size, r_z_array)

    if (small_int /= r_small_int .or. big_int /= r_big_int .or. small_real /= r_small_real .or. big_real /= r_big_real .or. z_size /= r_z_size) then
        write (*,*) 'rank ', rank, ' read did not match write'
    else
        write (*,*) 'rank ', rank, ' read matched write'
    endif

    call MPI_Barrier (MPI_COMM_WORLD, ierr)

    call adios_finalize (rank, ierr)

    call MPI_Finalize (ierr)
end program adios_test

subroutine test_write (group, filename, group_comm, small_int, big_int, small_real, big_real, a_size, a_array)
    implicit none
    include 'mpif.h'
    character (*), intent(in) :: group
    character (*), intent(in) :: filename
    integer, intent (in) :: group_comm
    integer, intent(in) :: small_int
    integer*8, intent(in) :: big_int
    real, intent(in) :: small_real
    real*8, intent(in) :: big_real
    integer, intent(in) :: a_size
    real, intent(in) :: a_array (a_size)

    integer :: istep1
    integer :: istep2
    integer :: istep3

    integer*8 :: handle
    integer*8 :: total_size
    integer :: err
    integer*8 :: size

    size = 900 * 1024

    istep1 = 11
    istep2 = 22
    istep3 = 33

    size = 4 + 4 + 8 + 4 + 8 + 4 + a_size * 4 + 4 + 4 + 4

    call adios_open (handle, trim(group)//char(0), trim(filename)//char(0), "w"//char(0), err)

    call adios_group_size (handle, size, total_size, group_comm, err)
    call adios_write (handle, "group_comm"//char(0), group_comm, err)

    call adios_write (handle, "small_int"//char(0), small_int, err)
    call adios_write (handle, "big_int"//char(0), big_int, err)
    call adios_write (handle, "small_real"//char(0), small_real, err)
    call adios_write (handle, "big_real"//char(0), big_real, err)
    call adios_write (handle, "ze0size"//char(0), a_size, err)
    call adios_write (handle, "zelectron0"//char(0), a_array, err)

    call adios_write (handle, "istep1"//char(0), istep1, err)
    call adios_write (handle, "istep2"//char(0), istep2, err)
    call adios_write (handle, "istep3"//char(0), istep3, err)

    call adios_close (handle, err)

end subroutine test_write

subroutine test_read (group, filename, group_comm, small_int, big_int, small_real, big_real, a_size, a_array)
    implicit none
    include 'mpif.h'
    character (*), intent(in) :: group
    character (*), intent(in) :: filename
    integer, intent (in) :: group_comm
    integer, intent(out) :: small_int
    integer*8, intent(out) :: big_int
    real, intent(out) :: small_real
    real*8, intent(out) :: big_real
    integer, intent(inout) :: a_size
    real, intent(out) :: a_array (a_size)

    integer*8 :: total_size
    integer :: err

    integer :: istep1
    integer :: istep2
    integer :: istep3

    integer*8 :: handle

    istep1 = 11
    istep2 = 22
    istep3 = 33

    call adios_open (handle, trim(group)//char(0), trim(filename)//char(0), "r"//char(0), err)

    call adios_group_size (handle, 0, total_size, 0, group_comm, err)
    call adios_read (handle, "small_int"//char(0), small_int, 4, err)
    call adios_read (handle, "big_int"//char(0), big_int, 8, err)
    call adios_read (handle, "small_real"//char(0), small_real, 4, err)
    call adios_read (handle, "big_real"//char(0), big_real, 8, err)
    call adios_read (handle, "ze0size"//char(0), a_size, 4, err)
    call adios_read (handle, "zelectron0"//char(0), a_array, a_size * 4, err)

    call adios_close (handle, err)

end subroutine test_read

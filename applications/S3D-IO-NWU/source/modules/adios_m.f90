!
! Author: Wei-keng Liao
!         wkliao@ece.northwestern.edu
!

!==============================================================================
module adios_m
    ! module for Read-Write Restart files using ADIOS

    implicit none
    include 'mpif.h'

    public :: adios_start, adios_wr, adios_end

    contains

!----< start_adios() >---------------------------------------------------
subroutine start_adios
    implicit none
    integer ierr
#ifdef _BUILD_ADIOS_F
    call adios_init ("config.xml"//char(0), ierr)
#endif
end subroutine start_adios

!----< end_adios() >-----------------------------------------------------
subroutine end_adios
    use topology_m, only : myid

    implicit none
    integer ierr
#ifdef _BUILD_ADIOS_F
    call adios_finalize (myid, ierr)
#endif
end subroutine end_adios

!----< write_adios() >------------------------------------------------------
subroutine write_adios(filename)
    use topology_m,  only : gcomm
    use param_m,     only : nsc, nx, ny, nz
    use variables_m, only : temp, pressure, yspecies, u
    use runtime_m,   only : time, tstep, time_save
    use bc_m,        only : pout
    use chemkin_m,   only : species_name
    use reference_m

    implicit none

    ! declarations passed in
    character*(*), intent(in) :: filename

    ! local variables
    integer i, ierr
    integer*8 :: handle
    integer*8 :: size, total_size
    character (len=200) :: group

#ifdef _BUILD_ADIOS_F
    group = "restart"

    call adios_open (handle, trim(group)//char(0), trim(filename)//char(0), "w"//char(0), ierr)

    size = nx * ny * nz * 8
    size = size * ((nsc+1) + 5) + 3 * 8 + 4 * 8

    call adios_group_size (handle, size, total_size, gcomm, ierr)

    ! local subarray dimensionality
    call adios_write (handle, "nx"//char(0), nx, ierr)
    call adios_write (handle, "ny"//char(0), ny, ierr)
    call adios_write (handle, "nz"//char(0), nz, ierr)

    call adios_write (handle,      "time"//char(0),      time, ierr)
    call adios_write (handle,     "tstep"//char(0),     tstep, ierr)
    call adios_write (handle, "time_save"//char(0), time_save, ierr)
    call adios_write (handle,      "pout"//char(0),      pout, ierr)

    !---- write array yspecies
    do i=1, nsc+1
        call adios_write (handle, "yspecies_"//trim(species_name(i))//char(0), yspecies(:,:,:,i), ierr)
    enddo

    !---- write array temp
    call adios_write (handle, "temp"//char(0), temp, ierr)

    !---- write array pressure
    call adios_write (handle, "pressure"//char(0), pressure, ierr)

    !---- write array u
    call adios_write (handle, "velocity_u"//char(0), u(:,:,:,1), ierr)
    call adios_write (handle, "velocity_v"//char(0), u(:,:,:,2), ierr)
    call adios_write (handle, "velocity_w"//char(0), u(:,:,:,3), ierr)

    call adios_close (handle, ierr)
#endif

    return
end subroutine write_adios

!==============================================================================
end module adios_m


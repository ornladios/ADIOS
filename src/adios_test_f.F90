!example program for using the Asynchronous IO system
program adios_test
    implicit none
    include 'mpif.h'
    integer*8 group1_id
    integer*8 io_handle1
    integer ierr
    integer ptracer_dim(1)
    integer partd_comm

    integer adios_integer, adios_real
    parameter (adios_integer=1, adios_real=2)

    integer istep 
    integer mype, nhybrid, istep, me, mi, ntracer, mflux
    real etracer, ptracer(4)
    integer mpsi
    real rdtemi(2), rdteme(2), pfluxpsi(2)
    integer phisize
    real phi(2), phip00(2), zonali(2), zonale(2)
    integer zionsize
    real zion(2), zion0(2)
    integer phisavesize
    real phisave(2)
    integer zelectronsize
    real zelectron(2)
    integer ze0size
    real zelectron0(2)
    integer nparam, mimax, mgrid, mzeta, memax

    ! The grouping node identifiers (use the MPI rank)
    ! (See the last subroutine in setup.F90)
    integer previous, current, next

    ptracer_dim(1) = 4
    istep=4
    mype = 1
    nhybrid = 6
    istep = 1
    me = 2
    mi = 2
    ntracer = 2
    mflux = 2
    etracer = 3
    ptracer(1) = 1.0
    ptracer(2) = 2.0
    ptracer(3) = 3.0
    ptracer(4) = 4.0
    mpsi = 2
    rdtemi(1) = 5.0
    rdtemi(2) = 6.0
    rdteme(1) = 7.0
    rdteme(2) = 8.0
    pfluxpsi(1) = 9.0
    pfluxpsi(2) = 10.0
    phisize = 2
    phi(1) = 11.0
    phi(2) = 12.0
    phip00(1) = 13.0
    phip00(2) = 14.0
    zonali(1) = 15.0
    zonali(2) = 16.0
    zonale(1) = 17.0
    zonale(2) = 18.0
    zionsize = 2
    zion(1) = 19.0
    zion(2) = 20.0
    zion0(1) = 21.0
    zion0(2) = 22.0
    phisavesize = 2
    phisave(1) = 23.0
    phisave(2) = 24.0
    zelectronsize = 2
    zelectron(1) = 25.0
    zelectron(2) = 26.0
    ze0size = 2
    zelectron0(1) = 27.0
    zelectron0(2) = 28.0
    nparam = 6
    mimax = 7
    mgrid = 8
    mzeta = 2
    memax = 10

    call MPI_Init (ierr)

    ! This is setup in GTC to be the communication within the group
    call MPI_Comm_dup (MPI_COMM_WORLD, partd_comm, ierr)
    previous = -1
    current = 0
    next = -1
    
    ! how to get the error messages out
    ! integer l, r, e
    ! character*1000 s
    ! call MPI_ERROR_STRING (ierr, s, r, e)

    ! setup the buffering once in startup
    ! (buffer size [MB], MPI pieces...) TEMP: Add the MPI_COMM_WORLD handle since C version incompatible
    call adios_init ("config_fortran.xml"//char(0), MPI_COMM_WORLD, MPI_COMM_SELF, MPI_INFO_NULL)

    ! get the group for use in writing and other operations (once at startup or any other time)
    call adios_get_group (group1_id, 'restart'//char(0))

    ! do our normal work for an interation

    ! open this nodes connection to the downstream consumer/storage of the data
    ! (handle, group id, filename, filepath)
    call adios_open (io_handle1, group1_id, 'restart.0'//char(0))

    ! write the group communicator
    call adios_write (io_handle1, "group_comm"//char(0), partd_comm)

    ! write a restart's data
    ! (handle, variable name, var)
    call adios_write (io_handle1, "nparam"//char(0), nparam)
    call adios_write (io_handle1, "mimax"//char(0), mimax)
    call adios_write (io_handle1, "mgrid"//char(0), mgrid)
    call adios_write (io_handle1, "mzeta"//char(0), mzeta)
    call adios_write (io_handle1, "memax"//char(0), memax)
    call adios_write (io_handle1, "mype"//char(0), mype)
    call adios_write (io_handle1, "nhybrid"//char(0), nhybrid)
    call adios_write (io_handle1, "istep"//char(0), istep)
    call adios_write (io_handle1, "me"//char(0), me)
    call adios_write (io_handle1, "mi"//char(0), mi)
    call adios_write (io_handle1, "ntracer"//char(0), ntracer)
    call adios_write (io_handle1, "mflux"//char(0), mflux)

    if (mype == 0) then
        call adios_write (io_handle1, "etracer"//char(0), etracer)
        call adios_write (io_handle1, "ptracer"//char(0), ptracer)
    endif

    call adios_write (io_handle1, "mpsi"//char(0), mpsi)
    call adios_write (io_handle1, "phisize"//char(0), phisize)
    call adios_write (io_handle1, "zionsize"//char(0), zionsize)
    call adios_write (io_handle1, "rdtemi"//char(0), rdtemi)
    call adios_write (io_handle1, "rdteme"//char(0), rdteme)
    call adios_write (io_handle1, "pfluxpsi"//char(0), pfluxpsi)
    call adios_write (io_handle1, "phi"//char(0), phi)
    call adios_write (io_handle1, "phip00"//char(0), phip00)
    call adios_write (io_handle1, "zonali"//char(0), zonali)
    call adios_write (io_handle1, "zonale"//char(0), zonale)
    call adios_write (io_handle1, "zion"//char(0), zion)
    call adios_write (io_handle1, "zion0"//char(0), zion0)

    if (nhybrid > 0) then
        call adios_write (io_handle1, "phisavesize"//char(0), phisavesize)
        call adios_write (io_handle1, "zelectronsize"//char(0), zelectronsize)
        call adios_write (io_handle1, "ze0size"//char(0), ze0size)
        call adios_write (io_handle1, "phisave"//char(0), phisave)
        call adios_write (io_handle1, "zelectron"//char(0), zelectron)
        call adios_write (io_handle1, "zelectron0"//char(0), zelectron0)
    endif

    ! make sure that anything that needs to be done
    ! at the end of a data write is completed.
    ! (handle)
    call adios_close (io_handle1)

    ! mark the end of an iteration for transmission timing
    call adios_end_iteration ()

    ! do work until another restart is written
    ! ...
    ! hint to data transfer that it is safe to communicate now
    call adios_start_calculation ()
    ! ...

    ! hint to data transfer that it is safe to communicate now
    call adios_stop_calculation ()

    ! once the simulation is done, make sure we have our IO finished
    call adios_finalize ()

    call MPI_Finalize (ierr)
end program adios_test

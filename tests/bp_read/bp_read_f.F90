program read_bp_f
    implicit none
    include "mpif.h"

    integer :: gcnt, vcnt, acnt, tstart, tstop, vrank, vtype, vtimed, ierr
    integer :: comm,i,j,k,l,m
    integer*8 :: fh, gh
    integer, dimension(10) :: dims, start, readsize
    integer, dimension(1000) :: var, vartrue
    character (len=100), dimension(50) :: vnamelist
    character (len=100), dimension(50) :: gnamelist
    
    call MPI_Init (ierr)
    comm = MPI_COMM_WORLD

    call adios_fopen (fh, "testbp_c.bp"//char(0), comm, ierr)
    call adios_inq_file (fh,gcnt,vcnt,acnt,tstart,tstop,gnamelist,ierr) 

    call adios_gopen (fh, gh, gnamelist(1), ierr) 
    call adios_inq_group(gh, vcnt, vnamelist, ierr)

    do i=1,vcnt 
        write (0,'(i3, ") ", a)') i, vnamelist(i)
    enddo
    ! vtype 
    ! 0 byte 
    ! 1 short 
    ! 2 integer
    ! 4 long 
    ! 5 real 
    ! 6 double
    ! 7 long double
    ! 9 string
    ! 10 complex
    ! 11 double_complex
    ! 50 unsigned_byte 
    ! 51 unsigned_short
    ! 52 unsigned_integer
    ! 54 unsigned_long
#if 0
    write (*,*) "name    ",  "  ndim    ", "    timed"
    call adios_inq_var (gh, vnamelist(1),vtype, vrank, vtimed, dims, ierr)
    write (*,*)"-----------------------------"
    write (*,"(a10, i3, i3)") vnamelist(1),vrank, vtimed
    call adios_inq_var (gh, vnamelist(2),vtype, vrank, vtimed, dims, ierr)
    write (*,*)"-----------------------------"
    write (*,"(a10, i3, i3)") vnamelist(2),vrank, vtimed
    call adios_inq_var (gh, vnamelist(3),vtype, vrank, vtimed, dims, ierr)
    write (*,*)"-----------------------------"
    write (*,"(a10, i3, i3)") vnamelist(3),vrank, vtimed
    start(1:10)=0
    readsize(1:10)=1 
    ! -------------- read 1D data ------------------
    call adios_inq_var (gh, "int_1D"//char(0), vtype, vrank, vtimed, dims, ierr)
    readsize(1)=dims(1)
    call adios_get_var (gh, "int_1D"//char(0), var, start, readsize, 1, ierr)
    write (*,*)"-----------------------------"
    write (*,"(a6, i3)") "int_1D", dims(1)
    write (*, "(10i3)") var(1:readsize(1))
    ! -------------- read 2D data ------------------
    call adios_inq_var (gh, "int_2D"//char(0),vtype, vrank, vtimed, dims, ierr)
    call adios_get_var (gh, "int_2D"//char(0), var, start, readsize, 1, ierr)
    write (*,*)"-----------------------------"
    write (*,"(a6, a3, 2i2, a3, 2i4, a1)") "int_2D","[",dims(1:2),"] [",readsize(1:2), "]"
    write (*,"(10i3)") var(1:readsize(1)*readsize(2))
    ! -------------- read 3D data ------------------
    call adios_inq_var (gh,"int_3D"//char(0),vtype, vrank, vtimed, dims, ierr)
    start(1:10)=0
    readsize(1)=dims(1) 
    readsize(2:10)=1
    readsize(2)=1
    call adios_get_var (gh, "int_3D"//char(0), var, start, readsize, 1, ierr)
    write (*,*)"-----------------------------"
    write (*,"(a6, a3, 3i3, a3, 3i3, a1)") "int_3D",'[',dims(1:3),'] [',readsize(1:3),"]"
    write (*,"(20i4)") var(1:readsize(1)*readsize(2)*readsize(3))
#endif

    ! -------------- read 4D data ------------------
    call adios_inq_var (gh, "int_4D"//char(0),vtype, vrank, vtimed, dims, ierr)
    readsize(1)=dims(1) 
    readsize(2)=1
    readsize(3)=1
    readsize(4)=2
    call adios_get_var (gh, "int_4D"//char(0), var, start, readsize, 1, ierr)
    write (*,*)"-----------------------------"
    write (*,"(a6, a2, 4i3, a3, 4i3, a1)") "int_4D","[",dims(1:4),"] [",readsize(1:4),"]" 
    write (*,"(80i5)") var(1:readsize(1)*readsize(2)*readsize(3)*readsize(4))
    vcnt = 1
    do i=0,readsize(1)-1
    do j=0,readsize(2)-1
    do k=0,readsize(3)-1
    do l=0,readsize(4)-1
        acnt = ((i*dims(2)+j)*dims(3)+k)*dims(4)+l
        vartrue(vcnt)=acnt
        vcnt = vcnt + 1
    enddo
    enddo
    enddo
    enddo
    write (*,"(80i5)") vartrue(1:readsize(1)*readsize(2)*readsize(3)*readsize(4))

#if 1
    ! -------------- read 5D data ------------------
    call adios_inq_var (gh, "int_5D"//char(0),vtype, vrank, vtimed, dims, ierr)
    readsize(1)=dims(1)
    readsize(2:10)=1
    readsize(3)=dims(3)

    vcnt = 1
    do i=0,readsize(1)-1
    do j=0,readsize(2)-1
    do k=0,readsize(3)-1
    do l=0,readsize(4)-1
    do m=0,readsize(5)-1
        acnt = (((i*dims(2)+j)*dims(3)+k)*dims(4)+l)*dims(5)+m
        vartrue(vcnt)=acnt
        vcnt = vcnt + 1
    enddo
    enddo
    enddo
    enddo
    enddo
    write (*,*)"-----------------------------"
    write (*,"(a, a2, 5i3, a3, 5i3, a1)") "int_5D","[",dims(1:5),"] [",readsize(1:5),"]" 
    call adios_get_var (gh, "int_5D"//char(0), var, start, readsize, 1, ierr)
    write (*,"(100i5)")var(1:readsize(1)*readsize(2)*readsize(3)*readsize(4)*readsize(5))
    write (*,"(100i5)")vartrue(1:readsize(1)*readsize(2)*readsize(3)*readsize(4)*readsize(5))
#endif
    call adios_gclose(gh, ierr)
    call adios_fclose(fh, ierr)

    call MPI_Finalize (ierr)
end program

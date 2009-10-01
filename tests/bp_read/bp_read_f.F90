program read_bp_f
    implicit none
    include "mpif.h"

    integer :: gcnt, vcnt, acnt, tstart, ntsteps, tstop, vrank, vtype, timedim, ierr
    integer :: comm,i,j,k,l,m
    integer*8 :: fh, gh
    integer*8, dimension(10) :: dims, start, readsize
    integer, dimension(1000) :: var, vartrue
    !character, dimension(:), allocatable :: varchar
    !integer, dimension(:), allocatable :: varint
    !real*8, dimension(:), allocatable :: vardouble
    character(100000) :: varchar
    integer, dimension(100000) :: varint
    real*8, dimension(10000000) :: vardouble
    integer :: totalsize
    character (len=100), dimension(5000) :: vnamelist
    character (len=100), dimension(5000) :: anamelist
    character (len=100), dimension(5000) :: gnamelist

    integer :: ilom,ihip,jlom,jhip,klom,khip
    integer, dimension(3) :: vstart, readcount
    integer,dimension(48) :: bconds  ! all boundary conditions as one array
    character(20)    :: vname
    integer          :: one
    real*8,dimension(102,66,3)  :: b1

    call MPI_Init (ierr)
    comm = MPI_COMM_WORLD

    varchar = ' '
    !call adios_fopen (fh, "/lustre/spider/scratch/pnorbert/TRACKP.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "TRACKP_00010.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "g_1x4_5x1.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "pixie3d.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "g_2x2_2x2_t1.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "xgcp.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "pgood.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "xgc.flowdiag.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "xgc.restart.000.03600.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "g.bp"//char(0), comm, ierr)
    !call adios_fopen (fh, "testbp.bp", comm, ierr)
    call adios_fopen (fh, "record.bp", comm, ierr)

    call adios_inq_file (fh,gcnt,vcnt,acnt,tstart,ntsteps,gnamelist,ierr) 
    tstop = ntsteps+ntsteps-1
    write (0,'("Number of groups : ",i0)') gcnt

    do i=1,gcnt 
        write (0,"(i5, a, a, a)") i,")  [", trim(gnamelist(i)),"] "
        write (10,"(i5, a, a, a)") i,")  [", trim(gnamelist(i)),"] "
    enddo

    call adios_gopen (fh, gh, gnamelist(1), ierr) 
    call adios_inq_group(gh, vcnt, vnamelist, acnt, anamelist, ierr)

    write (0,'("Number of variables in group ",a,": ",i0)') trim(gnamelist(1)), vcnt
    do i=1,vcnt 
        write (0,"(i5, a, a, a, i0)") i,")  [", trim(vnamelist(i)),"] ",index(vnamelist(i),"/var/")
        write (10,"(i5, a, a, a, i0)") i,")  [", trim(vnamelist(i)),"] ",index(vnamelist(i),"/var/")
    enddo

    write (0,'("Number of attributes in group ",a,": ",i0)') trim(gnamelist(1)), acnt
    do i=1,acnt 
        write (0,"(i5, a, a, a, i0)") i,")  [", trim(anamelist(i)),"] ",index(vnamelist(i),"/var/")
        write (10,"(i5, a, a, a, i0)") i,")  [", trim(anamelist(i)),"] ",index(vnamelist(i),"/var/")
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
    !write (*,*) "name    ",  "  ndim    ", "    dims"

    !start(1:10)=0
    !readsize(1:10)=1 
    !do i=tstart,tstop
    !    write(*,'("step ",i0,": gh = ",i0)') i, gh
    !    call adios_get_var (gh, "time"//char(0), varint, start, readsize, i, ierr)
    !    write(*,'("time[",i0,"] = ",i0)') i, varint(1)
    !enddo

    write (*,*)"-----------------------------"
    do i=1,vcnt 
        call adios_inq_var (gh, vnamelist(i), vtype, vrank, dims, timedim, ierr)
        start(1:10)=0
        readsize(1:10)=1 
        totalsize=1
        do j=1,vrank
            readsize(j) = dims(j)
            totalsize=totalsize*dims(j)
        enddo
        write (0,'(a," vrank=",i0," type=",i0," dim1=",i0," size=",i0)') trim(vnamelist(i)),vrank,vtype,dims(1),totalsize
        if (vtype == 0) then
            !if(allocated(varchar)) deallocate(varchar)
            !allocate(varchar(totalsize))
            !varchar(1)="!"
            call adios_read_var (gh, vnamelist(i), start, readsize, varchar, ierr)
        else if (vtype == 2) then
            !if(allocated(varint)) deallocate(varint)
            !allocate(varint(totalsize))
            !varint(1)=5
            !print *, varint(1)
            call adios_read_var (gh, vnamelist(i), start, readsize, varint, ierr)
            !print *, varint(1)
        else if (vtype == 6) then
            !if(allocated(vardouble)) deallocate(vardouble)
            !allocate(vardouble(totalsize))
            call adios_read_var (gh, vnamelist(i), start, readsize, vardouble, ierr)
        else
            write (0,'(a16,": Only integer or double type is handled here")') trim(vnamelist(i))
        endif
        if (vrank == 0) then
            if (vtype == 0) then
                write(*,'(a," = ",a)') trim(vnamelist(i)),varchar(1:1)
            else if (vtype == 2) then
                write(*,'(a," = ",i0)') trim(vnamelist(i)),varint(1)
                !print *, trim(vnamelist(i))
            else if (vtype == 6) then 
                write(*,'(a," = ",d20.10)') trim(vnamelist(i)),vardouble(1)
            endif
        else
            !write (*,*)"-----------------------------"
            write (*,'(a30,t32," [",i0,$)') vnamelist(i), dims(1)
            do j=2,vrank
                write (*,'(",",i0,$)') dims(j)
            enddo
            if (vtype == 0) then
                write (*,'("] = ",a)') trim(varchar)
            else
                write (*,'("]")') 
            endif
        endif
    enddo

#if 0
    call adios_get_var (gh, "bconds"//char(0), bconds, 0, 48, 1, ierr)
    write (*,'("ierr=",i0)') ierr
    write (*,'("bconds=",48i2)') bconds

    call adios_inq_var (gh, "/name/v1"//char(0), vtype, vrank, vtimed, dims, ierr)
    write (*,'("/name/v1: type=",i0," ndims=",i0)') vtype, vrank
    call adios_get_var (gh, "/name/v1"//char(0), vname, 0, 1, 1, ierr)
    write (*,'("ierr=",i0)') ierr
    write (*,'("/name/v1=[",a,"]")') trim(vname)

    ilom = 0
    jlom = 0
    klom = 0
    ihip = 101
    jhip = 65
    khip = 2
    vstart = (/ 0,0,0 /)
    !readcount = (/ 102,66,3 /)
    readcount = (/  ihip-ilom+1,jhip-jlom+1,khip-klom+1 /)
    call adios_get_var (gh, "/var/v1"//char(0), b1, vstart, readcount, 3, ierr)
    write (*,'("ierr=",i0)') ierr
    write (*,'("b1(1:10,0,0)=",10d20.12)') b1(1:10,1:1,1:1)
#endif

    call adios_gclose(gh, ierr)
    call adios_fclose(fh, ierr)

    call MPI_Finalize (ierr)
end program

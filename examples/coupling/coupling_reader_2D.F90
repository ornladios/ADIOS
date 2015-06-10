!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!
!  Coupling writer/reader
!
!  Read a 2D array 
!  to/from a file or from/to another coupling code using DataSpaces/ADIOS
!
!  npx * npy     processes write a 2D array, where each process writes an
!  
!  Data written
!    xy          2D global array as received from the writers
!    xy2         2D array with block,* decomp
!
!    All data are real*8 type
!
! (c) Oak Ridge National Laboratory, 2009
! Author: Norbert Podhorszki
!
module coupling_reader_2D_comm
    ! arguments
    integer :: npx, npy, npz  ! # of processors in x-y-z direction
    integer :: read_method    ! 0=bp, 1=hdf5, 2=dart, 3=dimes
    integer :: read_mode      ! 0=whole array on each process
                              ! 1=1D decomposition on 1st dim

    integer :: nproc          ! number of processors

    real*8, dimension(:,:),   allocatable :: xy, xy2

    ! Offsets and sizes
    integer :: ox, oy
    integer :: ldx, ldy
    integer :: gdx, gdy

    ! MPI variables
    integer :: group_comm, self_comm
    integer :: rank
    integer :: ierr

    ! ADIOS variables
    character (len=200) :: group
    !character (len=6)   :: nprocstr
    integer*8 :: adios_totalsize, adios_groupsize, adios_buf_size
    integer   :: adios_err
    integer*8 :: inh ! input file handle
    integer*8 :: adios_handle ! output file handle (used by gwrite too)
    character(len=256)  :: infn ! filename
    
    integer   :: wts  ! writer's output timestep index (read from 1,2...)
    logical, parameter   :: dump_text = .true.

end module coupling_reader_2D_comm


program coupling
    use coupling_reader_2D_comm
    use adios_read_mod
    implicit none
    include 'mpif.h'
    integer :: key, color
    real*8 :: t1,t2,io_time


    call MPI_Init (ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, key, ierr)

    ! Split MPI_COMM_WORLD for MPMD execution
    color = 2
    call MPI_Comm_split(MPI_COMM_WORLD, color, key, group_comm, ierr)
    call MPI_Comm_dup (MPI_COMM_SELF, self_comm, ierr)
    call MPI_Comm_rank (group_comm, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    call processArgs()
    if (rank == 0) then
        print '(" Process number        : ",i0," x ",i0)',  npx,npy
        print '(" Method for reading: ",i0)',  read_method

        if (nproc .ne. npx*npy) then
            print '(" Error: Number of processors ",i0,"does not match N*M=",i0)', nproc, npx*npy
            call exit(1)
        endif
    endif

    write (infn,'("writer.bp")')
    if (rank==0)    print *, " Input file: "//trim(infn)


    call adios_read_init_method (read_method, group_comm, "verbose=4; app_id=2; poll_interval=1000", ierr);
    call adios_init("coupling_reader_2D.xml", group_comm, ierr)
    !call MPI_Barrier (group_comm, ierr)

    !write (*,*) 'rank ', rank, "init completed"
    !write (nprocstr,'(i0)') nproc

    call adios_read_open (inh, infn, read_method, group_comm, ADIOS_LOCKMODE_CURRENT, 180.0, ierr)
    if (ierr .ne. 0) then
        print '(" Failed to open stream: ",a)', infn
        print '(" open stream ierr=: ",i0)', ierr
        call exit(1)
    endif

    t1=0;
    t2=0;
    call MPI_Barrier (group_comm, ierr)
    wts = 1  ! start reading from writer's 1st step
    do while (ierr==0)
        call MPI_Barrier (group_comm, ierr)
        t1=MPI_Wtime();
        call readArrays()
        call MPI_Barrier (group_comm, ierr)
        t2=t2+(MPI_Wtime()-t1);
        call adios_release_step(inh, ierr)
        call printArrays()
        call advanceArrays()
        call MPI_Barrier (group_comm, ierr)
        call adios_advance_step(inh, 0, 10.0, ierr)
        if (ierr==err_end_of_stream .and. rank==0) then
            print *, " Stream has terminated. Quit reader"
        elseif (ierr==err_step_notready .and. rank==0) then
            print *, " Next step has not arrived for a while. Assume termination" 
        endif
        wts = wts+1
    enddo
    call adios_read_close (inh, ierr)

!    call mpi_allreduce(t2, io_time, 1, MPI_INTEGER, MPI_MAX,group_comm, ierr)
    if(rank==0) then
        print '("Total read time = ", d12.2)', t2 
    endif

    ! Terminate
    call MPI_Barrier (group_comm, ierr)
    call adios_read_finalize_method (read_method, ierr)
    call adios_finalize (rank, ierr)
    call MPI_Finalize (ierr)
end program coupling



!!***************************
subroutine readArrays()
    use coupling_reader_2D_comm
    use adios_read_mod
    implicit none
    integer             :: i,j,k
    integer*8           :: sel
    integer*8, dimension(3) :: offset, readsize


    ! read gdx and gdy
    print '("rank=",i0,": Read in gdx and gdy, step",i0," from ",a)', rank, wts, trim(infn)
    call adios_get_scalar (inh,"gdx",gdx, ierr)
    call adios_get_scalar (inh,"gdy",gdy, ierr)
    print '("rank=",i0,": Got scalars gdx ",i0," gdy ",i0)', rank, gdx, gdy

    ! Calculate the local x,y offsets
    if (read_mode == 0) then
        ! each process reads the whole array
        ldx = gdx
        ldy = gdy
        ox = 0
        oy = 0
    elseif (read_mode == 1) then
        ldx = gdx / nproc 
        ldy = gdy
        ox = rank * ldx
        oy = 0
        ! last process should read all the rest of columns
        if (rank == nproc-1) then
            ldx = gdx - ldx*(nproc-1)
        endif
    endif


    if (.not. allocated(xy)) then
        allocate( xy  (1:ldx, 1:ldy) )
        allocate( xy2 (1:ldx, 1:ldy) )
    endif
    xy = -1.0 ! should be overwritten at reading the array next


    ! Read in the whole xy array
    readsize(1) = ldx
    readsize(2) = ldy
    offset(1)   = ox
    offset(2)   = oy
    call adios_selection_boundingbox (sel, 2, offset, readsize)
    print '("rank=",i0,": Read in xy(",i0,":",i0,",",i0,":",i0,") from ",a)', rank, &
            offset(1)+1, offset(1)+readsize(1), &
            offset(2)+1, offset(2)+readsize(2),trim(infn)
    call adios_schedule_read (inh, sel, "xy", 0, 1, xy, ierr)
    call adios_perform_reads (inh, ierr)
    call adios_selection_delete (sel)

    do j=1,ldy
        do i=1,ldx
            xy2(i,j) = rank*1.0 + xy(i,j)/100.0
        enddo
    enddo
    

    call MPI_Barrier (group_comm, ierr)

end subroutine readArrays

!!***************************
subroutine printArrays()
    use coupling_reader_2D_comm
    use adios_write_mod
    implicit none
    integer, parameter  :: u=20
    character(len=256)  :: outfn
    integer             :: i,j,k,writer

    if (read_mode == 0) then

        writer = mod(wts,nproc)
        if (writer == rank) then
    
            write (outfn, '("reader_",i3.3,".bp")') wts
            call adios_open (adios_handle, "reader2D", outfn, "w", self_comm, adios_err)
#include "gwrite_reader2D.fh"
            ! start streaming from buffer to disk
            call adios_close (adios_handle, adios_err)
            print '("rank=",i0,": One process wrote xy and new xy2 to ",a)', rank, trim(outfn)
    
        endif
    
        if (dump_text) then
    
            write (outfn, '("reader_",i3.3,"_",i0,".txt")') wts, rank
            open (u, FILE=outfn, STATUS='NEW', FORM="FORMATTED")
    
            ! print xyz_bbb
            write (u,'("xy(1:",i0,",1:",i0,") = ")') ldx, ldy
            do j=1,ldy
                do i=1,ldx
                    write (u, '(f6.2," ",$)') xy(i,j)
                enddo
                write (u,*) " "  ! new line
            enddo
    
            close (u)
    
        endif

    elseif (read_mode ==1) then
            write (outfn, '("reader_",i3.3,".bp")') wts
            call adios_open (adios_handle, "reader2D", outfn, "w", group_comm, adios_err)
#include "gwrite_reader2D.fh"
            ! start streaming from buffer to disk
            call adios_close (adios_handle, adios_err)
            print '("rank=",i0,": Collectively wrote xy and new xy2 to ",a)', rank, trim(outfn)
    endif

end subroutine printArrays

!!***************************
subroutine advanceArrays()
    use coupling_reader_2D_comm
    implicit none

    xy = 0.9 * xy
    xy2 = 0.9 * xy2

end subroutine advanceArrays


!!***************************
subroutine usage()
    print *, "Usage: coupling_reader_2D N M method read-mode"
    print *, "N:           number of processes in X dimension"
    print *, "M:           number of processes in Y dimension"
    print *, "read method: DATASPACES or DIMES for memory-to-memory coupling,"
    print *, "             otherwise file-based"
    print *, "read-mode:   0: each process reads whole global array"
    print *, "             1: 1D decomposition in 1st dimension"
end subroutine usage

!!***************************
subroutine processArgs()
    use coupling_reader_2D_comm
    use adios_read_mod

#ifndef __GFORTRAN__
!#ifndef __GNUC__
    interface
         integer function iargc()
         end function iargc
    end interface
!#endif
#endif

    character(len=256) :: npx_str, npy_str, method_str, mode_str, ts_str
    integer :: numargs

    !! process arguments
    numargs = iargc()
    !print *,"Number of arguments:",numargs
    if ( numargs < 4 ) then
        call usage()
        call exit(1)
    endif
    call getarg(1, npx_str)
    call getarg(2, npy_str)
    read (npx_str,'(i5)') npx
    read (npy_str,'(i5)') npy

    call getarg(3, method_str)
    if (trim(method_str) .eq. "DATASPACES") then
        read_method = ADIOS_READ_METHOD_DATASPACES 
    elseif (trim(method_str) .eq. "DIMES") then
        read_method = ADIOS_READ_METHOD_DIMES 
    elseif (trim(method_str) .eq. "FLEXPATH") then
        read_method = ADIOS_READ_METHOD_FLEXPATH 
    else
        read_method = ADIOS_READ_METHOD_BP
    endif

    call getarg(4, mode_str)
    read (mode_str,'(i5)') read_mode
    if (read_mode < 0 .or. read_mode > 1) then
        print *, "Argument error: read_mode must be 0 or 1"
        call usage()
        call exit(2)
    endif


end subroutine processArgs

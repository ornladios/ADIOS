!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!
!  GENARRAY
!
!  Write an ADIOS BP file from many processor for test purposes. 
!
!  nx * ny * nz     processes write a 3D array, where each process writes an
!  ndx * ndy * ndz  piece with filling with its rank as integer (4 bytes) value
!  
!
! (c) Oak Ridge National Laboratory, 2009
! Author: Norbert Podhorszki
!
module genarray_comm
    ! arguments
    character(len=256) :: outputfile, inputfile
    integer :: npx, npy, npz  ! # of processors in x-y-z direction
    integer :: ndx, ndy, ndz  ! size of array per processor
    logical :: common_size    ! .true.  if common local sizes are given as argument
                              ! .false. if we have to read sizes from a file

    integer :: gndx, gndy, gndz  ! size of the global array
    integer :: offx,offy,offz    ! offsets of local array in the global array

    integer, dimension(:,:,:), allocatable :: int_xyz

    ! MPI variables
    integer :: group_comm
    integer :: rank, nproc
    integer :: ierr

    ! ADIOS variables
    character (len=200) :: group
    character (len=200) :: filename
    !character (len=6)   :: nprocstr
    integer*8 :: handle, total_size, group_size
    integer   :: err

end module genarray_comm


program genarray
    use genarray_comm
    implicit none
    include 'mpif.h'

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    call adios_init ("genarray.xml", ierr)
    !call MPI_Barrier (group_comm, ierr)

    call processArgs()
    if (rank == 0) then
        print *,"Output file: "//trim(outputfile)
        print '(" Process number        : ",i0," x ",i0," x ",i0)', npx,npy,npz
        if (common_size) then
            print '(" Array size per process: ",i0," x ",i0," x ",i0)', ndx,ndy,ndz
        else
            print *," Array sizes per processes taken from file: "//trim(inputfile)
        endif

        if (nproc .ne. npx*npy*npz) then
            print '(" Error: Number of processors ",i0,"does not match ndx*ndy*ndz=",i0)', nproc, npx*npy*npz
            call exit(1)
        endif
    endif

    !write (*,*) 'rank ', rank, "init completed"
    !write (nprocstr,'(i0)') nproc

    call determineLocalSize()
    call determineGlobalSize()
    call determineOffsets()
    call generateLocalArray()
    call writeArray()

    ! Terminate
    call MPI_Barrier (MPI_COMM_WORLD, ierr)
    call adios_finalize (rank, ierr)
    call MPI_Finalize (ierr)
end program genarray


!!***************************
subroutine determineLocalSize()
    use genarray_comm
    implicit none
    if (common_size) then
       ! we are done since we know them from argument
    else
       ! have to read from file
       print *, "To be implemented: read sizes from file 1"
       call exit(2)
    endif
end subroutine determineLocalSize

!!***************************
subroutine determineGlobalSize()
    use genarray_comm
    implicit none
    if (common_size) then
        gndx = npx * ndx
        gndy = npy * ndy
        gndz = npz * ndz
    else
       ! have to read from file
       print *, "To be implemented: read sizes from file 2"
       call exit(2)
    endif
end subroutine determineGlobalSize

!!***************************
subroutine determineOffsets()
    use genarray_comm
    implicit none
    integer :: posx, posy, posz ! position index in the array
    if (common_size) then
        posx = mod(rank, npx)     ! 1st dim easy: 0, npx, 2npx... are in the same X position
        posy = mod(rank/npx, npy) ! 2nd dim: (0, npx-1) have the same dim (so divide with npx first)
        posz = rank/(npx*npy)     ! 3rd dim: npx*npy processes belong into one dim
        offx = posx * ndx
        offy = posy * ndy
        offz = posz * ndz
        print '("rank=",i0," pos: ",i0,",",i0,",",i0," offset: ",i0,",",i0,","i0)',  &
                rank, posx, posy, posz, offx, offy, offz
    else
       ! have to read from file
       print *, "To be implemented: read sizes from file 3"
       call exit(2)
    endif
end subroutine determineOffsets


!!***************************
subroutine generateLocalArray()
    use genarray_comm
    implicit none
    integer :: i,j,k
    allocate( int_xyz(1:ndx, 1:ndy, 1:ndz) )
    do k=1,ndz
        do j=1,ndy
            do i=1,ndx
                int_xyz(i,j,k) = rank 
            enddo
        enddo
    enddo
end subroutine generateLocalArray


!!***************************
subroutine writeArray()
    use genarray_comm
    implicit none
    ! Write out data using ADIOS
    group = "genarray"
    !  calculate how much we write from this processor
    group_size = 4 * 13              + &  ! X,Y,Z, nproc, all size_ and offs_ integers
                 4 * ndx * ndy * ndz      ! int_xyz

    !print '("rank=",i0," group=",A," file=",A," group_size=",i0)', rank, trim(group), &
    !    trim(outputfile), group_size
    call adios_open (handle, group, outputfile, "w", group_comm, err)
    call adios_group_size (handle, group_size, total_size, err)
    !print '("rank=",i0," total_size=",i0," err=",i0)', rank, total_size, err

    ! write dimensions and nproc 
    call adios_write (handle, "X", gndx, err)
    call adios_write (handle, "Y", gndy, err)
    call adios_write (handle, "Z", gndz, err)
    call adios_write (handle, "npx", npx, err)
    call adios_write (handle, "npy", npy, err)
    call adios_write (handle, "npz", npz, err)
    call adios_write (handle, "nproc", nproc, err)

    call adios_write (handle, "size_x", ndx, err)
    call adios_write (handle, "size_y", ndy, err)
    call adios_write (handle, "size_z", ndz, err)
    call adios_write (handle, "offs_x", offx, err) 
    call adios_write (handle, "offs_y", offy, err)
    call adios_write (handle, "offs_z", offz, err)
    call adios_write (handle, "int_xyz", int_xyz, err) 

    ! start streaming from buffer to disk
    call adios_close (handle, err)
    print '("rank=",i0,": write completed")', rank
end subroutine writeArray


!!***************************
subroutine usage()
    print *, "Usage: genarray  output N  M  K  [nx  ny  nz | infile]"
    print *, "output: name of output file"
    print *, "N:      number of processes in X dimension"
    print *, "M:      number of processes in Y dimension"
    print *, "K:      number of processes in Z dimension"
    print *, "nx:     local array size in X dimension per processor"
    print *, "ny:     local array size in Y dimension per processor"
    print *, "nz:     local array size in Z dimension per processor"
    print *, "infile: file that describes nx ny nz for each processor"
end subroutine usage

!!***************************
subroutine processArgs()
    use genarray_comm

#ifndef __GFORTRAN__
#ifndef __GNUC__
    interface
         integer function iargc()
         end function iargc
    end interface
#endif
#endif

    character(len=256) :: npx_str, npy_str, npz_str, ndx_str, ndy_str, ndz_str 
    integer :: numargs

    !! process arguments
    numargs = iargc()
    !print *,"Number of arguments:",numargs
    if ( numargs < 5 ) then
        call usage()
        call exit(1)
    endif
    call getarg(1, outputfile)
    call getarg(2, npx_str)
    call getarg(3, npy_str)
    call getarg(4, npz_str)
    read (npx_str,'(i5)') npx
    read (npy_str,'(i5)') npy
    read (npz_str,'(i5)') npz
    if ( numargs == 5 ) then
        call getarg(5, inputfile)
        ndx = 0
        ndy = 0
        ndz = 0
        common_size = .false.
    else if (numargs == 7) then
        call getarg(5, ndx_str)
        call getarg(6, ndy_str)
        call getarg(7, ndz_str)
        read (ndx_str,'(i6)') ndx
        read (ndy_str,'(i6)') ndy
        read (ndz_str,'(i6)') ndz
        inputfile=char(0)
        common_size = .true.
    else
        call usage()
        call exit(1)
    endif

end subroutine processArgs

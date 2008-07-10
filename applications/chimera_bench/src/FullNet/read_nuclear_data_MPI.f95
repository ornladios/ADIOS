SUBROUTINE read_nuclear_data_MPI( data_dir, data_desc )
!-----------------------------------------------------------------------------  
!  This routine reads, from the file netwinv, the nuclei included along 
!  with the nuclear data which will be needed for later calculations.  This 
!  data includes the atomic number, the number of protons and neutrons, and 
!  the binding energy (calculated from the tabulated mass excess).  Also the 
!  tabulations of the  partition functions, g, are read in for later 
!  interpolation. Once the set of nuclear data is read in, it is assigned 
!  to the proper nuclei.
!-----------------------------------------------------------------------------  

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE edit_module, ONLY : nlog
USE nuclear_data, ONLY : nname, aa, zz, nn, be
USE nuc_number, ONLY : ny
USE part_funct_data, ONLY : g, gg, t9i, angm
USE parallel_module, ONLY : myid, ierr

USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER (LEN=*),  INTENT(in)               :: data_dir

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

CHARACTER (LEN=80), INTENT(out)              :: data_desc

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                           :: var_name
CHARACTER (len=5)                            :: nam                                 
CHARACTER (len=5), ALLOCATABLE, DIMENSION(:) :: c_data    ! character array of initial data

INTEGER                                      :: i_extent  ! broadcast array extent
INTEGER, ALLOCATABLE, DIMENSION(:)           :: i_data    ! integer array of initial data
INTEGER                                      :: i,l,n,m,na,nb                                     
INTEGER                                      :: istat     ! allocation status

REAL(KIND=double)                            :: mex,a,sp
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: d_data

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in read_nuclear_data_MPI')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_nuclear_data_MPI')

!-----------------------------------------------------------------------
!  Allocate buffer
!-----------------------------------------------------------------------

ALLOCATE (i_data(1), STAT = istat)
IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  character (LEN=80) :: data_dir,data_desc
!  data_dir='../../FullNet/Data_alpha'
!
!  Read in the data description
!-----------------------------------------------------------------------

  OPEN (13,file=TRIM(data_dir)//"/net_desc",STATUS='old')
  READ (13,"(A)") data_desc
  OPEN (12,FILE=TRIM(data_dir)//"/netwinv",STATUS='old')
  READ (12,"(i5)") ny

!-----------------------------------------------------------------------
!  Load buffer
!-----------------------------------------------------------------------

  i_data(1)              = ny

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast array dimension to all processors, deallocate buffer.
!-----------------------------------------------------------------------

i_extent                 = 1
CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( i_data, i_extent, MPI_INTEGER           , 0, MPI_COMM_WORLD, ierr)
ny                       = i_data(1)

DEALLOCATE (i_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate data buffers
!-----------------------------------------------------------------------

ALLOCATE (c_data(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c_data    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (i_data(24), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_data(29*ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,1001) var_name; END IF

i_data                   = 0
d_data                   = zero

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Read in the partition function iteration grid, and fix endpoints
!-----------------------------------------------------------------------

  READ (12,"(24i3)") (i_data(i),i=1,24)

  DO  n= 1,ny
    READ (12,"(a5)") nam
  END DO

!-----------------------------------------------------------------------
!  Read in nuclear parameters and partition function interpolation table
!   into buffers for broadcast to all processors.
!-----------------------------------------------------------------------

  DO l = 1,ny
    READ (12,"(a5,f12.3,2i4,f6.1,f10.3)") nam,a,na,nb,sp,mex
    READ (12,"(8f9.2)") (d_data(29*(l-1)+m),m=1,24)
    d_data(29*(l-1)+25)  = a
    d_data(29*(l-1)+26)  = dble(float(na))
    d_data(29*(l-1)+27)  = dble(float(nb))
    d_data(29*(l-1)+28)  = 8.07144 * dble(float(nb)) + 7.28899 * dble(float(na)) - mex
    d_data(29*(l-1)+29)  = 2.d0 * sp + 1.d0
    c_data(l)            = nam
  END DO

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

i_extent                 = 5 * ny
CALL MPI_BCAST( c_data, i_extent, MPI_CHARACTER         , 0, MPI_COMM_WORLD, ierr)

i_extent                 = 24
CALL MPI_BCAST( i_data, i_extent, MPI_INTEGER           , 0, MPI_COMM_WORLD, ierr)

i_extent                 = 29 * ny
CALL MPI_BCAST( d_data, i_extent, MPI_DOUBLE_PRECISION  , 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Set size of nuclear parameter arrays and allocate.
!-----------------------------------------------------------------------

ALLOCATE (nname(0:ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nname     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aa(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aa        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zz(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zz        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nn(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nn        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (g(24,ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'g         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gg(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gg        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (angm(0:ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'angm      '; WRITE (nlog,1001) var_name; END IF

aa            = zero
zz            = zero
nn            = zero
be            = zero
g             = zero
gg            = zero
angm          = zero

!-----------------------------------------------------------------------
!  nname(0), gg(0) and angm(0) are placeholders for non-nuclei.
!-----------------------------------------------------------------------

nname(0)                 = '   '
angm(0)                  = 0.0d0

!-----------------------------------------------------------------------
!  Transfer nuclear parameters and partition function interpolation
!   table.
!-----------------------------------------------------------------------

DO i = 1,24
  t9i(i)                 = i_data(i) * 0.01d0
END DO
t9i(24)                  = t9i(24) * 10.d0

DO l = 1,ny
  DO m = 1,24
    g(m,l)               = d_data(29*(l-1)+m)
  END DO
  aa(l)                  = d_data(29*(l-1)+25)
  zz(l)                  = d_data(29*(l-1)+26)
  nn(l)                  = d_data(29*(l-1)+27)
  be(l)                  = d_data(29*(l-1)+28)
  angm(l)                = d_data(29*(l-1)+29)
  nname(l)               = c_data(l)
END DO

!-----------------------------------------------------------------------
!  Deallocate arrays.
!-----------------------------------------------------------------------

DEALLOCATE (c_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c_data    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (i_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (d_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE read_nuclear_data_MPI

SUBROUTINE read_reaction_data_MPI( data_dir )
!-----------------------------------------------------------------------------  
!  This routine reads in the necessary reaction data
!-----------------------------------------------------------------------------  

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE edit_module, ONLY : nlog
USE nuc_number, ONLY : ny
USE nuclear_data, ONLY : nname
USE ffn_data, ONLY : nffn, rf, r1, r2, ffnsum, ffnenu
USE cross_sect_data, ONLY : nreac, csect1, rc1, q1, rpf1, n1i, iwk1, ires1, &
& irev1, iffn, csect2, rc2, q2, rpf2, h2, n2i, iwk2, ires2, irev2, csect3, &
& rc3, q3, n3i, iwk3, ires3, irev3
USE parallel_module, ONLY : myid, ierr
USE reac_rate_data, ONLY : la, le, nan, mu1, a1, b1, n11, mu2, a2, b2, n21, &
& n22, mu3, a3, b3, n31, n32, n33

USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER (LEN=*),  INTENT(in)                 :: data_dir

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                             :: var_name
CHARACTER (len=1)                              :: null
CHARACTER (len=5), ALLOCATABLE, DIMENSION(:)   :: c_data   ! character array of initial data

INTEGER                                        :: i_extent  ! broadcast array extent
INTEGER                                        :: istat     ! allocation status
INTEGER                                        :: i,j,n,l
INTEGER                                        :: nr,idel,jdel,ijd,ndi,nds
INTEGER, ALLOCATABLE, DIMENSION(:)             :: i_data    ! integer array of initial data

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: ffnsum_t
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: ffnenu_t
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: d_data

      common /ceign/ nds, ndi,idel,jdel,ijd

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in read_nuclear_data_MPI')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_nuclear_data_MPI')

!-----------------------------------------------------------------------
!  Allocate buffer
!-----------------------------------------------------------------------

ALLOCATE (i_data(5), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,1001) var_name; END IF

i_data                           = 0

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  data_dir='../../FullNet/Data_alpha'
!  Read in nuclear set and numbers of reactions
!-----------------------------------------------------------------------

  OPEN (4,file=TRIM(data_dir)//"/nets4",FORM='unformatted',STATUS='old')
  READ(4) ny                                                               
  READ(4) (nname(i),i=1,ny)
  READ(4) nffn 
  READ(4) (nreac(i),i=1,3)

!-----------------------------------------------------------------------
!  Load buffer
!-----------------------------------------------------------------------

  i_data(1)                      = ny
  i_data(2)                      = nffn
  i_data(3:5)                    = nreac(1:3)

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

i_extent                        = 5
CALL MPI_BCAST( I_data, i_extent, MPI_INTEGER           , 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Transfer data to module arrays and deallocate buffer.
!-----------------------------------------------------------------------

ny                               = i_data(1)
nffn                             = i_data(2)
nreac(1:3)                       = i_data(3:5)

DEALLOCATE (i_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate buffer arrays
!-----------------------------------------------------------------------

ALLOCATE (c_data(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c_data    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Load buffer
!-----------------------------------------------------------------------

  c_data(1:ny)                   = nname(1:ny)

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

i_extent                         = 5 * ny
CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( c_data, i_extent, MPI_CHARACTER         , 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Transfer data to module arrays and deaalocate buffer.
!-----------------------------------------------------------------------

nname(1:ny)                      = c_data(1:ny)

DEALLOCATE (c_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c_data    '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  If there are FFN rates, read in the FFN data and set FFN array sizes
!-----------------------------------------------------------------------

IF ( nffn > 0 ) THEN

!-----------------------------------------------------------------------
!  Allocate buffer array
!-----------------------------------------------------------------------

  ALLOCATE (d_data(2*nffn*143), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,1001) var_name; END IF

  d_data                         = zero

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

  IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Read in the FFN data and set FFN array sizes
!-----------------------------------------------------------------------

    ALLOCATE (ffnsum_t(nffn,143), STAT = istat)
      IF ( istat /= 0 ) THEN; var_name = 'ffnsum_t  '; WRITE (nlog,1001) var_name; END IF
    ALLOCATE (ffnenu_t(nffn,143), STAT = istat)
      IF ( istat /= 0 ) THEN; var_name = 'ffnenu_t  '; WRITE (nlog,1001) var_name; END IF

    ffnsum_t                     = zero
    ffnenu_t                     = zero

    OPEN(3,FILE=TRIM(data_dir)//"/netweak",STATUS='old')
    DO i=1,nffn
      READ(3,"(a1)") null  
      READ(3,"(9(f8.3))") (ffnsum_t(i,j),ffnenu_t(i,j),j=1,143)  
    END DO

!-----------------------------------------------------------------------
!  Load buffers and deallocate temporary data arrays
!-----------------------------------------------------------------------

    DO i=1,nffn
      DO j = 1,143
        d_data(143*(i-1)+j)          = ffnsum_t(i,j)
        d_data(143*nffn+143*(i-1)+j) = ffnenu_t(i,j)
      END DO ! j
    END DO ! i

    DEALLOCATE (ffnsum_t, STAT = istat)
      IF ( istat /= 0 ) THEN; var_name = 'ffnsum_t  '; WRITE (nlog,2001) var_name; END IF
    DEALLOCATE (ffnenu_t, STAT = istat)
      IF ( istat /= 0 ) THEN; var_name = 'ffnenu_t  '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

  END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

  i_extent                       = 2 * nffn * 143
  CALL MPI_BCAST( d_data, i_extent, MPI_DOUBLE_PRECISION  , 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Allocate module arrays, transfer data, and deaallocate data buffers
!-----------------------------------------------------------------------

  ALLOCATE (rf(nffn), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'rf        '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (r1(nffn), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'r1        '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (r2(nffn), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'r2        '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ffnsum(nffn,143), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'ffnsum    '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ffnenu(nffn,143), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'ffnenu    '; WRITE (nlog,1001) var_name; END IF

  rf                             = zero
  r1                             = zero
  r2                             = zero
  ffnsum                         = zero
  ffnenu                         = zero

  DO i=1,nffn
    DO j = 1,143
      ffnsum(i,j)                = d_data(143*(i-1)+j)
      ffnenu(i,j)                = d_data(143*nffn+143*(i-1)+j)
    END DO ! j
  END DO ! i

  DEALLOCATE (d_data, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,2001) var_name; END IF

END IF ! nffn > 0

!-----------------------------------------------------------------------
!  Allocate buffer array for cross section data
!-----------------------------------------------------------------------

nr                               = nreac(1)
ALLOCATE (i_data(6*nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_data(8*nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,1001) var_name; END IF

i_data                           = 0
d_data                           = zero

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Allocate and read in reaction arrays for 1csect1 reactant reactions
!-----------------------------------------------------------------------

  OPEN(2,file=TRIM(data_dir)//"/nets3",FORM='unformatted',STATUS='old')

  ALLOCATE (rc1(7,nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'rc1       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (q1(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'q1        '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (n1i(3,nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'n1i       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (iwk1(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'iwk1      '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ires1(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'ires1     '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (irev1(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'irev1     '; WRITE (nlog,1001) var_name; END IF

  rc1                            = zero
  q1                             = zero
  n1i                            = 0
  iwk1                           = 0
  ires1                          = 0
  irev1                          = 0

  DO j = 1,nr
    READ(2) n,(n1i(l,j),l=1,3),iwk1(j),ires1(j),irev1(j),(rc1(l,j),l=1,7),q1(j)
    IF( n /= j ) THEN
      WRITE(6,*) 'Error in nets3, 1',j,n
      EXIT
    END IF
  END DO

!-----------------------------------------------------------------------
!  Load buffers and deallocate temporary data arrays
!-----------------------------------------------------------------------

  DO j = 1,nr
    i_data(j)                    = n1i(1,j)
    i_data(nr+j)                 = n1i(2,j)
    i_data(2*nr+j)               = n1i(3,j)
    i_data(3*nr+j)               = iwk1(j)
    i_data(4*nr+j)               = ires1(j)
    i_data(5*nr+j)               = irev1(j)
    DO l = 1,7
      d_data(nr*(l-1)+j)         = rc1(l,j)
    END DO ! l
    d_data(nr*7+j)               = q1(j)
  END DO ! j

  DEALLOCATE (rc1, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'rc1       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (q1, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'q1        '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (n1i, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'n1i       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (iwk1, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'iwk1      '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (ires1, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'ires1     '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (irev1, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'irev1     '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

i_extent                         = 6 * nr
CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( i_data, i_extent, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)

i_extent                         = 8 * nr
CALL MPI_BCAST( d_data, i_extent, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Allocate and load module arrays, and deallocate buffers
!-----------------------------------------------------------------------

ALLOCATE (csect1(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'csect1    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rc1(7,nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rc1       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (q1(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'q1        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rpf1(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpf1      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n1i(3,nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n1i       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iwk1(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iwk1      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ires1(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ires1     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (irev1(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'irev1     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iffn(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iffn      '; WRITE (nlog,1001) var_name; END IF

csect1                           = zero
rc1                              = zero
q1                               = zero
rpf1                             = zero
n1i                              = 0
iwk1                             = 0
ires1                            = 0
irev1                            = 0
iffn                             = 0

DO j = 1,nr
  n1i(1,j)                       = i_data(j)
  n1i(2,j)                       = i_data(nr+j)
  n1i(3,j)                       = i_data(2*nr+j)
  iwk1(j)                        = i_data(3*nr+j)
  ires1(j)                       = i_data(4*nr+j)
  irev1(j)                       = i_data(5*nr+j)
  DO l = 1,7
    rc1(l,j)                     = d_data(nr*(l-1)+j)
  END DO ! l
  q1(j)                          = d_data(nr*7+j)
END DO ! j

DEALLOCATE (i_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (d_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate buffer array for cross section data
!-----------------------------------------------------------------------

nr                               = nreac(2)
ALLOCATE (i_data(7*nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_data(8*nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,1001) var_name; END IF

i_data                           = 0
d_data                           = zero

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Allocate and read in reaction arrays for 2 reactant reactions
!-----------------------------------------------------------------------

  ALLOCATE (rc2(7,nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'rc2       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (q2(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'q2        '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (n2i(4,nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'n2i       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (iwk2(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'iwk2      '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ires2(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'ires2     '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (irev2(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'irev2     '; WRITE (nlog,1001) var_name; END IF

  rc2                            = zero
  q2                             = zero
  n2i                            = 0
  iwk2                           = 0
  ires2                          = 0
  irev2                          = 0

  DO j = 1,nr
    READ(2)  n,(n2i(l,j),l=1,4),iwk2(j),ires2(j),irev2(j),(rc2(l,j),l=1,7),q2(j)                           
    IF( n /= j ) THEN
      WRITE(6,*) 'Error in nets3, 2',j,n
      EXIT
    END IF ! n /= j
  END DO ! j

!-----------------------------------------------------------------------
!  Load buffers and deallocate temporary data arrays
!-----------------------------------------------------------------------

  DO j = 1,nr
    i_data(j)                    = n2i(1,j)
    i_data(nr+j)                 = n2i(2,j)
    i_data(2*nr+j)               = n2i(3,j)
    i_data(3*nr+j)               = n2i(4,j)
    i_data(4*nr+j)               = iwk2(j)
    i_data(5*nr+j)               = ires2(j)
    i_data(6*nr+j)               = irev2(j)
    DO l = 1,7
      d_data(nr*(l-1)+j)         = rc2(l,j)
    END DO ! l
    d_data(nr*7+j)               = q2(j)
  END DO ! j

  DEALLOCATE (rc2, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'rc2       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (q2, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'q2        '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (n2i, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'n2i       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (iwk2, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'iwk2      '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (ires2, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'ires2     '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (irev2, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'irev2     '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

i_extent                         = 7 * nr
CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( i_data, i_extent, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)

i_extent                         = 8 * nr
CALL MPI_BCAST( d_data, i_extent, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Allocate and load module arrays, and deallocate buffers
!-----------------------------------------------------------------------

ALLOCATE (csect2(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'csect2    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rc2(7,nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rc2       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (q2(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'q2        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rpf2(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpf2      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (h2(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'h2        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n2i(4,nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n2i       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iwk2(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iwk2      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ires2(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ires2     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (irev2(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'irev2     '; WRITE (nlog,1001) var_name; END IF

csect2                           = zero
rc2                              = zero
q2                               = zero
rpf2                             = zero
n2i                              = 0
iwk2                             = 0
ires2                            = 0
irev2                            = 0

DO j = 1,nr
  n2i(1,j)                       = i_data(j)
  n2i(2,j)                       = i_data(nr+j)
  n2i(3,j)                       = i_data(2*nr+j)
  n2i(4,j)                       = i_data(3*nr+j)
  iwk2(j)                        = i_data(4*nr+j)
  ires2(j)                       = i_data(5*nr+j)
  irev2(j)                       = i_data(6*nr+j)
  DO l = 1,7
    rc2(l,j)                     = d_data(nr*(l-1)+j)
  END DO ! l
  q2(j)                          = d_data(nr*7+j)
END DO ! j

DEALLOCATE (i_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (d_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate buffer array reaction arrays for 2 reactant reactions
!-----------------------------------------------------------------------

nr                               = nreac(3)
ALLOCATE (i_data(6*nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_data(8*nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,1001) var_name; END IF

i_data                           = 0
d_data                           = zero

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Allocate and read in reaction arrays for 3 reactant reactions
!-----------------------------------------------------------------------

  ALLOCATE (rc3(7,nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'rc3       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (q3(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'q3        '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (n3i(3,nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'n3i       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (iwk3(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'iwk3      '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ires3(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'ires3     '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (irev3(nr), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'irev3     '; WRITE (nlog,1001) var_name; END IF

  rc3                            = zero
  q3                             = zero
  n3i                            = 0
  iwk3                           = 0
  ires3                          = 0
  irev3                          = 0

  DO j = 1,nr
    READ(2)  n,(n3i(l,j),l=1,3),iwk3(j),ires3(j),irev3(j),(rc3(l,j),l=1,7),q3(j)                                                
    IF(n/=j) THEN
      WRITE(6,*) 'Error in nets3, 3',j,n
      EXIT
    END IF
  END DO

!-----------------------------------------------------------------------
!  Load buffers and deallocate temporary data arrays
!-----------------------------------------------------------------------

  DO j = 1,nr
    i_data(j)                    = n3i(1,j)
    i_data(nr+j)                 = n3i(2,j)
    i_data(2*nr+j)               = n3i(3,j)
    i_data(3*nr+j)               = iwk3(j)
    i_data(4*nr+j)               = ires3(j)
    i_data(5*nr+j)               = irev3(j)
    DO l = 1,7
      d_data(nr*(l-1)+j)         = rc3(l,j)
    END DO ! l
    d_data(nr*7+j)               = q3(j)
  END DO ! j
      
  DEALLOCATE (rc3, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'rc3       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (q3, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'q3        '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (n3i, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'n3i       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (iwk3, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'iwk3      '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (ires3, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'ires3     '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (irev3, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'irev3     '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

i_extent                         = 6 * nr
CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( i_data, i_extent, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)

i_extent                         = 8 * nr
CALL MPI_BCAST( d_data, i_extent, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Allocate and load module arrays, and deallocate buffers
!-----------------------------------------------------------------------

ALLOCATE (csect3(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'csect3    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rc3(7,nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rc3       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (q3(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'q3        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n3i(3,nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n3i       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iwk3(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iwk3      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ires3(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ires3     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (irev3(nr), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'irev3     '; WRITE (nlog,1001) var_name; END IF

csect3                           = zero
rc3                              = zero
q3                               = zero
n3i                              = 0
iwk3                             = 0
ires3                            = 0
irev3                            = 0

DO j = 1,nr
  n3i(1,j)                       = i_data(j)
  n3i(2,j)                       = i_data(nr+j)
  n3i(3,j)                       = i_data(2*nr+j)
  iwk3(j)                        = i_data(3*nr+j)
  ires3(j)                       = i_data(4*nr+j)
  irev3(j)                       = i_data(5*nr+j)
  DO l = 1,7
    rc3(l,j)                     = d_data(nr*(l-1)+j)
  END DO ! l
  q3(j)                          = d_data(nr*7+j)
END DO ! j

DEALLOCATE (i_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (d_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate module arrays and buffer
!-----------------------------------------------------------------------

ALLOCATE (la(3,ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'la        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (le(3,ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'le        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (i_data(6*ny+5), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,1001) var_name; END IF

la                               = 0
le                               = 0
i_data                           = 0

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Read in the data linking nuclei to the reactions which  affect them.
!   Also read in the matrix sparseness descriptors.
!-----------------------------------------------------------------------

  DO i = 1,ny
    READ(4)  n,la(1,i),le(1,i),la(2,i),le(2,i),la(3,i),le(3,i)
    IF( n /= i) THEN
      WRITE(6,*) 'Error in nets4',i,n
    END IF
  END DO
  READ(4) idel,jdel,ijd,ndi,nds

!-----------------------------------------------------------------------
!  Load buffer
!-----------------------------------------------------------------------

  DO i = 1,ny
    i_data(i)                    = la(1,i)
    i_data(ny+i)                 = la(2,i)
    i_data(2*ny+i)               = la(3,i)
    i_data(3*ny+i)               = le(1,i)
    i_data(4*ny+i)               = le(2,i)
    i_data(5*ny+i)               = le(3,i)
  END DO
  i_data(6*ny+1)                 = idel
  i_data(6*ny+2)                 = jdel
  i_data(6*ny+3)                 = ijd
  i_data(6*ny+4)                 = ndi
  i_data(6*ny+5)                 = nds

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

i_extent                         = 6 * ny + 5
CALL MPI_BCAST( i_data, i_extent, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Load module arrays, and deallocate buffer
!-----------------------------------------------------------------------

DO i = 1,ny
  la(1,i)                        = i_data(i)
  la(2,i)                        = i_data(ny+i)
  la(3,i)                        = i_data(2*ny+i)
  le(1,i)                        = i_data(3*ny+i) 
  le(2,i)                        = i_data(4*ny+i)
  le(3,i)                        = i_data(5*ny+i)
END DO
idel                             = i_data(6*ny+1)
jdel                             = i_data(6*ny+2)
ijd                              = i_data(6*ny+3)
ndi                              = i_data(6*ny+4)
nds                              = i_data(6*ny+5)

DEALLOCATE (i_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate module arrays and buffers
!-----------------------------------------------------------------------

nan(1)      =le(1,ny)
ALLOCATE (mu1(nan(1)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'mu1       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a1(nan(1)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a1        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (b1(nan(1)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b1        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n11(nan(1)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n11       '; WRITE (nlog,1001) var_name; END IF

mu1                              = 0
a1                               = zero
b1                               = zero
n11                              = 0

nan(2)      = le(2,ny)
ALLOCATE (mu2(nan(2)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'mu2       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a2(nan(2)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a2        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (b2(nan(2)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b2        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n21(nan(2)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n21       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n22(nan(2)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n22       '; WRITE (nlog,1001) var_name; END IF

mu2                              = 0
a2                               = zero
b2                               = zero
n21                              = 0
n22                              = 0

nan(3)      = le(3,ny)
ALLOCATE (mu3(nan(3)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'mu3       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a3(nan(3)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a3        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (b3(nan(3)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b3        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n31(nan(3)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n31       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n32(nan(3)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n32       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n33(nan(3)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n33       '; WRITE (nlog,1001) var_name; END IF

mu3                              = 0
a3                               = zero
b3                               = zero
n31                              = 0
n32                              = 0
n33                              = 0

ALLOCATE (i_data(nan(1)+nan(2)+nan(3)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (d_data(nan(1)+nan(2)+nan(3)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,1001) var_name; END IF

i_data                           = 0
d_data                           = zero

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Read extended reaction->nuclei arrays links
!-----------------------------------------------------------------------

  DO j = 1,nan(1)
    READ( 2) a1(j),mu1(j) 
  END DO

  DO j=1,nan(2)
    READ( 2)  a2(j),mu2(j)  
  END DO

  DO j = 1,nan(3)
    READ( 2)  a3(j),mu3(j)
  END DO

!-----------------------------------------------------------------------
!  Load buffers
!-----------------------------------------------------------------------

  DO j = 1,nan(1)
    i_data(j)                    = mu1(j)
    d_data(j)                    = a1(j)
  END DO

  DO j = 1,nan(2)
    i_data(nan(1)+j)             = mu2(j)
    d_data(nan(1)+j)             = a2(j)
  END DO

  DO j = 1,nan(3)
    i_data(nan(1)+nan(2)+j)      = mu3(j)
    d_data(nan(1)+nan(2)+j)      = a3(j)
  END DO

!-----------------------------------------------------------------------
!              ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Broadcast data to all processors.
!-----------------------------------------------------------------------

i_extent                         = nan(1) + nan(2) + nan(3)
CALL MPI_BCAST( i_data, i_extent, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)

i_extent                         = nan(1) + nan(2) + nan(3)
CALL MPI_BCAST( d_data, i_extent, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!-----------------------------------------------------------------------
!  Load module arrays, and deallocate buffers
!-----------------------------------------------------------------------

DO j = 1,nan(1)
  mu1(j)                         = i_data(j)
  a1(j)                          = d_data(j)
END DO

DO j = 1,nan(2)
  mu2(j)                         = i_data(nan(1)+j)
  a2(j)                          = d_data(nan(1)+j)
END DO

DO j = 1,nan(3)
  mu3(j)                         = i_data(nan(1)+nan(2)+j)
  a3(j)                          = d_data(nan(1)+nan(2)+j)
END DO

DEALLOCATE (i_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'i_data    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (d_data, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'd_data    '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Fill extended reaction->nuclei arrays
!-----------------------------------------------------------------------




DO j = 1,nan(1)
  n11(j)                         = n1i(1,mu1(j))
END DO

DO j=1,nan(2)
  n21(j)                         = n2i(1,mu2(j))
  n22(j)                         = n2i(2,mu2(j))
END DO

DO j = 1,nan(3)
  n31(j)                         = n3i(1,mu3(j))
  n32(j)                         = n3i(2,mu3(j))
  n33(j)                         = n3i(3,mu3(j))
END DO

RETURN
END SUBROUTINE read_reaction_data_MPI

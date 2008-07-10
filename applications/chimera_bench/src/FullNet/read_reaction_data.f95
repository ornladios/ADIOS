SUBROUTINE read_reaction_data( data_dir )
!-----------------------------------------------------------------------------  
!  This routine reads in the necessary reaction data
!-----------------------------------------------------------------------------  

USE numerical_module, ONLY: zero

USE nuc_number, ONLY : ny
USE nuclear_data, ONLY : nname
USE ffn_data, ONLY : nffn, rf, r1, r2, ffnsum, ffnenu
USE cross_sect_data, ONLY : nreac, csect1, rc1, q1, rpf1, n1i, iwk1, ires1, &
& irev1, iffn, csect2, rc2, q2, rpf2, h2, n2i, iwk2, ires2, irev2, csect3, &
& rc3, q3, n3i, iwk3, ires3, irev3
USE reac_rate_data, ONLY : la, le, nan, mu1, a1, b1, n11, mu2, a2, b2, n21, &
& n22, mu3, a3, b3, n31, n32, n33

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER (LEN=*),  INTENT(in)  :: data_dir

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=1) :: null

INTEGER       :: myid = 12
INTEGER       :: i,j,n,l
INTEGER       :: nr,idel,jdel,ijd,ndi,nds

      common /ceign/ nds, ndi,idel,jdel,ijd

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
!  IF there are FFN rates, read in the FFN data and set FFN array sizes
!-----------------------------------------------------------------------

IF ( nffn > 0 ) THEN
  OPEN(3,FILE=TRIM(data_dir)//"/netweak",STATUS='old')
  ALLOCATE (rf(nffn),r1(nffn),r2(nffn))
  ALLOCATE (ffnsum(nffn,143),ffnenu(nffn,143))

  rf        = zero
  r1        = zero
  r2        = zero
  ffnsum    = zero
  ffnenu    = zero

  DO i=1,nffn
    READ(3,"(a1)") null  
    READ(3,"(9(f8.3))") (ffnsum(i,j),ffnenu(i,j),j=1,143)  
  END DO
END IF

!-----------------------------------------------------------------------
!  Read in the reaction cross section data
!-----------------------------------------------------------------------

OPEN(2,file=TRIM(data_dir)//"/nets3",FORM='unformatted',STATUS='old')

!-----------------------------------------------------------------------
!  Allocate and read in reaction arrays for 1 reactant reactions
!-----------------------------------------------------------------------

nr          = nreac(1)
ALLOCATE (csect1(nr),rc1(7,nr),q1(nr),rpf1(nr))
ALLOCATE (n1i(3,nr),iwk1(nr),ires1(nr),irev1(nr),iffn(nr))

csect1      = zero
rc1         = zero
q1          = zero
rpf1        = zero
n1i         = 0
iwk1        = 0
ires1       = 0
irev1       = 0
iffn        = 0

DO j = 1,nr
  READ(2) n,(n1i(l,j),l=1,3),iwk1(j),ires1(j),irev1(j),(rc1(l,j),l=1,7),q1(j)  
  IF( n /= j ) THEN
    WRITE(6,*) 'Error in nets3, 1',j,n
    EXIT
  END IF
END DO

!-----------------------------------------------------------------------
!  Allocate and read in reaction arrays for 2 reactant reactions
!-----------------------------------------------------------------------

nr          = nreac(2)
ALLOCATE (csect2(nr),rc2(7,nr),q2(nr),rpf2(nr),h2(nr))
ALLOCATE (n2i(4,nr),iwk2(nr),ires2(nr),irev2(nr))

csect2      = zero
rc2         = zero
q2          = zero
rpf2        = zero
n2i         = 0
iwk2        = 0
ires2       = 0
irev2       = 0

DO j = 1,nr
  READ(2)  n,(n2i(l,j),l=1,4),iwk2(j),ires2(j),irev2(j),(rc2(l,j),l=1,7),q2(j)                           
  IF( n /= j ) THEN
    WRITE(6,*) 'Error in nets3, 2',j,n
    EXIT
  END IF
END DO

!-----------------------------------------------------------------------
!  Allocate and read in reaction arrays for 3 reactant reactions
!-----------------------------------------------------------------------

nr          = nreac(3)
ALLOCATE (csect3(nr),rc3(7,nr),q3(nr))
ALLOCATE (n3i(3,nr),iwk3(nr),ires3(nr),irev3(nr))

csect3      = zero
rc3         = zero
q3          = zero
n3i         = 0
iwk3        = 0
ires3       = 0
irev3       = 0

DO j = 1,nr
  READ(2)  n,(n3i(l,j),l=1,3),iwk3(j),ires3(j),irev3(j),(rc3(l,j),l=1,7),q3(j)                                                
  IF(n/=j) THEN
    WRITE(6,*) 'Error in nets3, 3',j,n
    EXIT
  END IF
END DO

!-----------------------------------------------------------------------
!  Allocate and read in the data linking nuclei to the reactions which 
!  affect them.  Also read in the matrix sparseness descriptors.
!-----------------------------------------------------------------------

ALLOCATE (la(3,ny),le(3,ny))

la          = 0
le          = 0

DO i = 1,ny                                                             

  READ(4)  n,la(1,i),le(1,i),la(2,i),le(2,i),la(3,i),le(3,i)
  IF( n /= i) THEN
    WRITE(6,*) 'Error in nets4',i,n
  END IF

END DO
READ(4) idel,jdel,ijd,ndi,nds

!-----------------------------------------------------------------------
!  Create and fill extended reaction->nuclei arrays
!-----------------------------------------------------------------------

nan(1)      = le(1,ny)
ALLOCATE (mu1(nan(1)),a1(nan(1)),b1(nan(1)),n11(nan(1)))

mu1         = 0
a1          = zero
b1          = zero
n11         = 0

DO j = 1,nan(1)
  READ( 2) a1(j),mu1(j) 
  n11(j)    = n1i(1,mu1(j))
END DO

nan(2)      = le(2,ny)
ALLOCATE (mu2(nan(2)),a2(nan(2)),b2(nan(2)))
ALLOCATE (n21(nan(2)),n22(nan(2)))

mu2         = 0
a2          = zero
b2          = zero
n21         = 0
n22         = 0

DO j=1,nan(2)
  READ( 2)  a2(j),mu2(j)  
  n21(j)    = n2i(1,mu2(j))
  n22(j)    = n2i(2,mu2(j))
END DO

nan(3)      = le(3,ny)
ALLOCATE (mu3(nan(3)),a3(nan(3)),b3(nan(3)))
ALLOCATE (n31(nan(3)),n32(nan(3)),n33(nan(3)))

mu3         = 0
a3          = zero
b3          = zero
n31         = 0
n32         = 0
n33         = 0

DO j = 1,nan(3)
  READ( 2)  a3(j),mu3(j)
  n31(j)    = n3i(1,mu3(j))
  n32(j)    = n3i(2,mu3(j))
  n33(j)    = n3i(3,mu3(j))
END DO

RETURN
END SUBROUTINE read_reaction_data

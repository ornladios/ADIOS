MODULE Ecuyer_random
! L'Ecuyer's 1996 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! N.B. This version is compatible with Lahey's ELF90
! http://www.ozemail.com.au/~milleraj
! Latest revision - 30 March 1999

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

! These are unsigned integers in the C version
INTEGER, SAVE :: s1 = 1234, s2 = -4567, s3 = 7890

CONTAINS

SUBROUTINE init_seeds(i1, i2, i3)
IMPLICIT NONE

INTEGER, INTENT(IN) :: i1, i2, i3

s1 = i1
s2 = i2
s3 = i3
IF (IAND(s1,-2) == 0) s1 = i1 - 1023
IF (IAND(s2,-8) == 0) s2 = i2 - 1023
IF (IAND(s3,-16) == 0) s3 = i3 - 1023

RETURN
END SUBROUTINE init_seeds



FUNCTION taus88() RESULT(random_numb)
! Generates a random number between 0 and 1.  Translated from C function in:
! Reference:
! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
! generators', Math. of Comput., 65, 203-213.

! The cycle length is claimed to be about 2^(88) or about 3E+26.
! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

IMPLICIT NONE
REAL (dp) :: random_numb

INTEGER   :: b

! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
!      to the left if j > 0, otherwise to the right.

b  = ISHFT( IEOR( ISHFT(s1,13), s1), -19)
s1 = IEOR( ISHFT( IAND(s1,-2), 12), b)
b  = ISHFT( IEOR( ISHFT(s2,2), s2), -25)
s2 = IEOR( ISHFT( IAND(s2,-8), 4), b)
b  = ISHFT( IEOR( ISHFT(s3,3), s3), -11)
s3 = IEOR( ISHFT( IAND(s3,-16), 17), b)
random_numb = IEOR( IEOR(s1,s2), s3) * 2.3283064365E-10_dp + 0.5_dp

RETURN
END FUNCTION taus88

END MODULE Ecuyer_random

C
C***********************************************************************
C
C    MODULE:       MATSCL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Multiply a N by M Matrix by a scalar
C
C    CALL LINE:    CALL MATSCL(A,SCALAR,B,N,M)
C
C    INPUTS:       A = Array to be scaled  (D)
C                  SCALAR = Constant to multiply matrix by (D)
C                  N = Number of rows in arrays (I)
C                  M = Number of columns in arrays (I)
C
C    OUTPUTS:      B = Array containing SCALAR x A (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATSCL(A,SCALAR,B,N,M)
C
      IMPLICIT NONE
C
      INTEGER N, M
      DOUBLE PRECISION A(N,M), B(N,M), SCALAR
C
C
C                 Local variables
C
      INTEGER I, J
C
C                 Loop over all elements of the array
      DO 20 J=1,M,1
        DO 10 I=1,N,1
C
          B(I,J) = SCALAR*A(I,J)
C
 10     CONTINUE
 20   CONTINUE
C
 999  RETURN
C
      END

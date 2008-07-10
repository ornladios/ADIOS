C
C***********************************************************************
C
C    MODULE:       MATMUL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Multiplies two Matrices (LxM)x(MxN)
C
C    CALL LINE:    CALL MATMUL(A,B,C,L,M,N)
C
C    INPUTS:       A,B= Arrays to be added  (D)
C                  L,M,N = Dimensions of arrays (I)
C
C    OUTPUTS:      C = Array containing A x B (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATMUL(A,B,C,L,M,N)
C
      IMPLICIT NONE
C
      INTEGER L, M, N
      DOUBLE PRECISION A(L,M), B(M,N), C(L,N)
C
C
C                 Local variables
C
      INTEGER I, J, K
      DOUBLE PRECISION SUM
C
C                 Loop over all elements of the array
      DO 30 I=1,L,1
        DO 20 J=1,N,1
C
C                 Initialize SUM for a new element
          SUM = 0.0
C                 Calculate (i,j)th element
          DO 10 K=1,M,1
            SUM = SUM+A(I,K)*B(K,J)
 10       CONTINUE
          C(I,J) = SUM
C
 20     CONTINUE
 30   CONTINUE
C
 999  RETURN
C
      END

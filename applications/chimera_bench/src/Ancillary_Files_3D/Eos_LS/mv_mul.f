C***********************************************************************
C
C    MODULE:       MV_MUL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Multiplies a Matrix times a vector (NxN)x(N)
C
C    CALL LINE:    CALL MV_MUL(A,V,RV,N)
C
C    INPUTS:       A = Array to be multiplied  (D)
C                  V = Vector to be multiplied (D)
C                  N = Dimensions of arrays & vector (I)
C
C    OUTPUTS:      RV = resultant vector (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MV_MUL(A,V,RV,N)
C
      IMPLICIT NONE
C
      INTEGER N
      DOUBLE PRECISION A(N,N), V(N), RV(N)
C
C
C                 Local variables
C
      INTEGER I, J
      DOUBLE PRECISION SUM
C
C                 Loop over all elements of the array
      DO 20 I=1,N,1
C
C                 Initialize SUM for a new element
        SUM = 0.0
C                 Calculate (i)th element
        DO 10 J=1,N,1
          SUM = SUM+A(I,J)*V(J)
 10     CONTINUE
        RV(I) = SUM
C
 20   CONTINUE
C
 999  RETURN
C
      END

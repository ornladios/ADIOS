C
C***********************************************************************
C
C    MODULE:       MATCOP
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Copy one N by M Matrix into another
C
C    CALL LINE:    CALL MATCOP(A,B,N,M)
C
C    INPUTS:       A = Array to be copied  (D)
C                  N = Number of rows in arrays (I)
C                  M = Number of columns in arrays (I)
C
C    OUTPUTS:      C = Array to be copied into (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATCOP(A,B,N,M)
C
      IMPLICIT NONE
C
      INTEGER N,M
      DOUBLE PRECISION A(N,M), B(N,M)
C
C
C                 Local variables
C
      INTEGER I, J
C
      DO 20 J=1,M,1
        DO 10 I=1,N,1
          B(I,J) = A(I,J)
 10     CONTINUE
 20   CONTINUE
C
 999  RETURN
C
      END

C***********************************************************************
C
C    MODULE:       MATADD
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Adds two N by M Matrices
C
C    CALL LINE:    CALL MATINV(A,B,C,N,M)
C
C    INPUTS:       A,B= Arrays to be added  (D)
C                  N = Number of rows in arrays (I)
C                  M = Number of columns in arrays (I)
C
C    OUTPUTS:      C = Array containing A+B (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATADD(A,B,C,N,M)
C
      IMPLICIT NONE
C
      INTEGER N,M
      DOUBLE PRECISION A(N,M), B(N,M), C(N,M)
C
C
C                 Local variables
C
      INTEGER I, J
C
      DO 20 J=1,M,1
        DO 10 I=1,N,1
          C(I,J) = A(I,J)+B(I,J)
 10     CONTINUE
 20   CONTINUE
C
 999  RETURN
C
      END

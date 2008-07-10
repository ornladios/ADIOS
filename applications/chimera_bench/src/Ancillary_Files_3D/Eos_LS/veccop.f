C
C
C***********************************************************************
C
C    MODULE:       VECCOP
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/10/90
C
C    PURPOSE:      Copy one vector on length N into another
C
C    CALL LINE:    CALL VECCOP(A,B,N)
C
C    INPUTS:       A = Array to be copied  (D)
C                  N = Number of rows in arrays (I)
C
C    OUTPUTS:      B = Array to be copied into (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE VECCOP(A,B,N)
C
      IMPLICIT NONE
C
      INTEGER N
      DOUBLE PRECISION A(N), B(N)
C
C
C                 Local variables
C
      INTEGER I
C
      DO 10 I=1,N,1
        B(I) = A(I)
 10   CONTINUE
C
 999  RETURN
C
      END

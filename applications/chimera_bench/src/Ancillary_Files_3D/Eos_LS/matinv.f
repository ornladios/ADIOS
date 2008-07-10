C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         DMATRIX.FOR
C
C***********************************************************************
C
C    MODULE:       MATINV
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Inverts a N by N Matrix
C
C    CALL LINE:    CALL MATINV(A,AINV,N)
C
C    INPUTS:       A = Array to be inverted  (D)
C                  N = dimesion of arrays (I)
C
C    OUTPUTS:      AINV = Inverse of A (D)
C
C    CALLS :       Numerical recipes routines LUDCMP, LUBKSB
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATINV(A,AINV,N)
C
      IMPLICIT NONE
C
      INTEGER N
      DOUBLE PRECISION A(N,N), AINV(N,N)
C
C
C                 Local variables
C
      INTEGER NPHYS, I, J
      PARAMETER(NPHYS=10)
      DOUBLE PRECISION TEMP(NPHYS,NPHYS), Y(NPHYS,NPHYS), D
      INTEGER INDEX(NPHYS)
C
C                 Make a copy of the array, and initialize
C                 the indentity matrix
      DO 20 J=1,N,1
        DO 10 I=1,N,1
          Y(I,J) = 0.0
          TEMP(I,J) = A(I,J)
 10     CONTINUE
        Y(J,J) = 1.0
 20   CONTINUE
C
C
C                 LU decompose the matrix
      CALL LUDCMP(TEMP,N,NPHYS,INDEX,D)
C
C
C                 Back substitute to get inverse
      DO 30 J=1,N,1
        CALL LUBKSB(TEMP,N,NPHYS,INDEX,Y(1,J))
 30   CONTINUE
C
C
C                 Copy temporary array into the inverse array
      DO 50 J=1,N,1
        DO 40 I=1,N,1
          AINV(I,J) = Y(I,J)
 40     CONTINUE
 50   CONTINUE
C
C
 999  RETURN
      END

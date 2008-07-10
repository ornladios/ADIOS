C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       F_3_2
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         11/29/89
C
C    CALL LINE:    F_3_2(Y)      (3/2th Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       3/2th Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION F_3_2(y)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(n=201)
      DIMENSION A(7),eta(n),f32(n),f32a(n),f12(n),fr(n),f12a(n),fra(n),
     > fia(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
      DATA A,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0
     1,6.41500299D-02,.4D0,1.32934039D0,1,n/
      IF(y .gt. 30.) goto 10
      if(y .lt. -10.) goto 20
      call splint(eta,f32,f32a,n,y,f1,klo,khi)
      GO TO 100
 10   X2=y**(-2)
      F1=A(6)*SQRT(y)*(1./X2+A(1)-(A(2)+X2*A(3))*X2)
      GO TO 100
 20   F0=DEXP(y)
      F1=A(7)*F0*(1.-(A(4)-(A(5)-.03125*F0)*F0)*F0)
 100  F_3_2=F1
 999  RETURN
      END

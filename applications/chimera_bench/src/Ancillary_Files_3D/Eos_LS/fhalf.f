C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       FHALF
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         3/10/90
C
C    CALL LINE:    FHALF(Y)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       Ratio of 1/2th Fermi Integral to the -1/2th Fermi
C                  Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION FHALF(y)
      IMPLICIT real*8(a-h,o-z)
      parameter(n=201)
      DIMENSION A(7),f12(n),fia(n),eta(n),f32(n),fr(n),f32a(n),
     > fra(n),f12a(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
      DATA A,th,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0
     1,6.41500299D-02,.4D0,1.32934039D0,.3333333333333d0,1,n/
      IF(y .gt. 30.) goto 10
      if(y .lt. -10.) goto 20
      call splint(eta,fr,fra,n,y,f1,klo,khi)
      GO TO 100
 10   X2=y**(-2)
      F1=y*th*(1.+(.2*A(1)+(.6*A(2)+1.4*X2*A(3))*X2)*x2)
     > /(1.-(.2*th*a(1)-(a(2)-4.2*x2*a(3))*x2)*x2)
      GO TO 100
 20   F0=EXP(y)
      F1=(1.-(2*a(4)-(3*a(5)-.125*f0)*f0)*f0)/
     > (2.-(8*a(4)-(18*a(5)-f0)*f0)*f0)
 100  FHALF=F1
 999  RETURN
      END

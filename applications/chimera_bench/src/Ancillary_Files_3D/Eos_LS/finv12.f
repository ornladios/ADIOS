C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       FINV12
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         11/29/89
C
C    CALL LINE:    FINV12(Y)      (Inverse of the 1/2th Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       Inverse of Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION FINV12(y)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(n=201)
      DIMENSION AI(8),f12(n),fia(n),eta(n),f32(n),fr(n),f32a(n),
     > fra(n),f12a(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
      DATA AI,klo,khi/-.822467032D0,-1.21761363D0,-9.16138616D0,
     1.398942281D0,.0732748216D0,-1.310707D0,1.12837917D0,
     28.2810645D-3,1,n/
      if(y .gt. 109.695) goto 10
      if(y .lt. 4.0234e-5) goto 20
      call splint(f12,eta,fia,n,y,f1,klo,khi)
      GO TO 100
 10   X2=(1.5*y)**(.666666667)
      X4=1./(X2*X2)
      F1=X2*(1.+(AI(1)+(AI(2)+AI(3)*X4)*X4)*X4)
      GO TO 100
 20   F1=LOG(AI(7)*MAX(y,1.D-20)*(1.+(AI(4)+(AI(5)+AI(8)*y)*y)*y))
 100  finv12=F1
 999  RETURN
      END

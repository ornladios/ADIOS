        SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,KLO,KHI)
c     Computes spline fit of Y(X); YA(XA) is input data, Y2A are spline
c  coefficents, klo and khi are running indices which bracket X.
        implicit real*8(a-h,o-z)
        DIMENSION XA(N),YA(N),Y2A(N)
cc  Determine the bracketing indices
        IF((KHI.GT.N).OR.(KLO.LT.1)) THEN
          WRITE(*,*) ' SPLINT CRASHED '
          WRITE(*,*) ' N,KHI,KLO =',N,KHI,KLO
          STOP
        ENDIF
C
        IF((KHI.GT.N).OR.(KHI.LT.1)) KHI = N
        IF((KLO.GT.N).OR.(KLO.LT.1)) KLO = 1
C
        IF(KHI-KLO.GT.1)GOTO 1
c-----------------------------------------------------
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2
        KHI=MAX(KHI-1,1)
        KLO=MAX(KLO-1,1)
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2
        KHI=MIN(KHI+2,N)
        KLO=MIN(KLO+2,N)
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2
c-----------------------------------------------------
        KLO=1
        KHI=N
 1      IF(KHI-KLO.EQ.1) GOTO 2
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
        KHI=K
        ELSE
        KLO=K
        ENDIF
        GOTO 1
 2      H=XA(KHI)-XA(KLO)
        IF(H.EQ.0.) PAUSE 'BAD XA INPUT. '
cc  Compute spline fit.
        A=(XA(KHI)-X)/H
        B=(X-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+
     1 ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*H**2/6.
c       write(*,5)klo,khi,x,xa(klo),xa(khi),ya(klo),ya(khi)
c     > ,y2a(klo),y2a(khi),y
 5      format(2i3,8(1pe9.2))
        RETURN
        END

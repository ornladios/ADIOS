        SUBROUTINE SPLINE(X,Y,N,Y2)
c  Computes spline coefficients; Y(X) is input data; Y2 is output.
        implicit real*8(a-h,o-z)
        DIMENSION X(N),Y(N),Y2(N),U(500)
        Y2(1)=0.
        U(1)=0.
        DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
 11     U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     1 /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
        Y2(N)=0.
        DO 12 K=N-1,1,-1
 12     Y2(K)=Y2(K)*Y2(K+1)+U(K)
        RETURN
        END

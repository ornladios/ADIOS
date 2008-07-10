      subroutine my_dcsint (ndata, xdata, fdata, break, cscoef)
      implicit none
      integer ndata, i, ibcbeg, ibcend
      real (kind=8) xdata(ndata), fdata(ndata), break(ndata), &
                    cscoef(4,ndata)
      do i=1,ndata
        break(i) = xdata(i)
        cscoef(1,i) = fdata(i)
      enddo
      ibcbeg = 0
      ibcend = 0
      call cubspl (break, cscoef, ndata, ibcbeg, ibcend)
      return
      end subroutine my_dcsint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function my_dppder(ideriv, x, korder, nintv, break, ppcoef)
      implicit none
      integer ideriv, korder, nintv
      real (kind=8) my_dppder, x, break(nintv+1), ppcoef(korder,nintv)
      real (kind=8) ppvalu
      my_dppder = ppvalu(break, ppcoef, nintv, korder, x, ideriv)
      return
      end function my_dppder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
!  from  * a practical guide to splines *  by c. de boor    
!     ************************  input  ***************************
!     n = number of data points. assumed to be .ge. 2.
!     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
!        data points. tau is assumed to be strictly increasing.
!     ibcbeg, ibcend = boundary condition indicators, and
!     c(2,1), c(2,n) = boundary condition information. specifically,
!        ibcbeg = 0  means no boundary condition at tau(1) is given.
!           in this case, the not-a-knot condition is used, i.e. the
!           jump in the third derivative across tau(2) is forced to
!           zero, thus the first and the second cubic polynomial pieces
!           are made to coincide.)
!        ibcbeg = 1  means that the slope at tau(1) is made to equal
!           c(2,1), supplied by input.
!        ibcbeg = 2  means that the second derivative at tau(1) is
!           made to equal c(2,1), supplied by input.
!        ibcend = 0, 1, or 2 has analogous meaning concerning the
!           boundary condition at tau(n), with the additional infor-
!           mation taken from c(2,n).
!     ***********************  output  **************************
!     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
!        of the cubic interpolating spline with interior knots (or
!        joints) tau(2), ..., tau(n-1). precisely, in the interval
!        (tau(i), tau(i+1)), the spline f is given by
!           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
!        where h = x - tau(i). the function program *ppvalu* may be
!        used to evaluate f or its derivatives from tau,c, l = n-1,
!        and k=4.
!
      implicit none
      integer ibcbeg,ibcend,n,   i,j,l,m
      real (kind=8) c(4,n),tau(n),   divdf1,divdf3,dtau,g
!
!****** a tridiagonal linear system for the unknown slopes s(i) of
!  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
!  ination, with s(i) ending up in c(2,i), all i.
!     c(3,.) and c(4,.) are used initially for temporary storage.
      l = n-1
!compute first differences of tau sequence and store in c(3,.). also,
!compute first divided difference of data and store in c(4,.).
      do 10 m=2,n
         c(3,m) = tau(m) - tau(m-1)
 10          c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
!construct first equation from the boundary condition, of the form
!             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
      if (ibcbeg-1)                     11,15,16
 11    if (n .gt. 2)                     go to 12
!     no condition at left end and n = 2.
      c(4,1) = 1.
      c(3,1) = 1.
      c(2,1) = 2.*c(4,2)
                                        go to 25
!     not-a-knot condition at left end and n .gt. 2.
 12                                      c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3,3)
      c(2,1) =((c(3,2)+2.*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3))/c(3,1)
                                        go to 19
!     slope prescribed at left end.
 15                                      c(4,1) = 1.
      c(3,1) = 0.
                                        go to 18
!     second derivative prescribed at left end.
 16                                      c(4,1) = 2.
      c(3,1) = 1.
      c(2,1) = 3.*c(4,2) - c(3,2)/2.*c(2,1)
 18    if(n .eq. 2)                      go to 25
!  if there are interior knots, generate the corresp. equations and car-
!  ry out the forward pass of gauss elimination, after which the m-th
!  equation reads    c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
 19     do 20 m=2,l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
 20          c(4,m) = g*c(3,m-1) + 2.*(c(3,m) + c(3,m+1))
!construct last equation from the second boundary condition, of the form
!           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
!     if slope is prescribed at right end, one can go directly to back-
!     substitution, since c array happens to be set up just right for it
!     at this point.
      if (ibcend-1)                     21,30,24
 21    if (n .eq. 3 .and. ibcbeg .eq. 0) go to 22
!     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
!     left end point.
      g = c(3,n-1) + c(3,n)
      c(2,n) = ((c(3,n)+2.*g)*c(4,n)*c(3,n-1) &
             + c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
      g = -g/c(4,n-1)
      c(4,n) = c(3,n-1)
                                        go to 29
!     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
!     knot at left end point).
 22                                      c(2,n) = 2.*c(4,n)
      c(4,n) = 1.
                                        go to 28
!     second derivative prescribed at right endpoint.
 24                                      c(2,n) = 3.*c(4,n) + c(3,n)/2.*c(2,n)
      c(4,n) = 2.
                                        go to 28
 25                                      if (ibcend-1)                     26,30,24
 26                                       if (ibcbeg .gt. 0)                go to 22
!     not-a-knot at right endpoint and at left endpoint and n = 2.
      c(2,n) = c(4,n)
                                        go to 30
 28                                      g = -1./c(4,n-1)
!complete forward pass of gauss elimination.
 29                                       c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
!carry out back substitution
 30    j = l 
 40        c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
         j = j - 1
         if (j .gt. 0)                  go to 40
!****** generate cubic coefficients in each interval, i.e., the deriv.s
!  at its left endpoint, from value and slope at its endpoints.
      do 50 i=2,n
         dtau = c(3,i)
         divdf1 = (c(1,i) - c(1,i-1))/dtau
         divdf3 = c(2,i-1) + c(2,i) - 2.*divdf1
         c(3,i-1) = 2.*(divdf1 - c(2,i-1) - divdf3)/dtau
 50          c(4,i-1) = (divdf3/dtau)*(6./dtau)
                                        return
      end subroutine cubspl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function ppvalu (break, coef, l, k, x, jderiv )
!  from  * a practical guide to splines *  by c. de boor    
!calls  interv
!calculates value at  x  of  jderiv-th derivative of pp fct from pp-repr
!
!******  i n p u t  ******
!  break, coef, l, k.....forms the pp-representation of the function  f
!        to be evaluated. specifically, the j-th derivative of  f  is
!        given by
!
!     (d**j)f(x) = coef(j+1,i) + h*(coef(j+2,i) + h*( ... (coef(k-1,i) +
!                             + h*coef(k,i)/(k-j-1))/(k-j-2) ... )/2)/1
!
!        with  h = x - break(i),  and
!
!       i  =  max( 1 , max( j ,  break(j) .le. x , 1 .le. j .le. l ) ).
!
!  x.....the point at which to evaluate.
!  jderiv.....integer giving the order of the derivative to be evaluat-
!        ed.  a s s u m e d  to be zero or positive.
!
!******  o u t p u t  ******
!  ppvalu.....the value of the (jderiv)-th derivative of  f  at  x.
!
!******  m e t h o d  ******
!     the interval index  i , appropriate for  x , is found through a
!  call to  interv . the formula above for the  jderiv-th derivative
!  of  f  is then evaluated (by nested multiplication).
!
      implicit none
      integer jderiv,k,l,   i,m,ndummy
      real (kind=8) ppvalu, break(l+1),coef(k,l),x,   fmmjdr,h
      ppvalu = 0.
      fmmjdr = k - jderiv
!              derivatives of order  k  or higher are identically zero.
      if (fmmjdr .le. 0.)               go to 99
!
!              find index  i  of largest breakpoint to the left of  x .
      call interv ( break, l+1, x, i, ndummy )
!
!      Evaluate  jderiv-th derivative of  i-th polynomial piece at  x .
      h = x - break(i)
      m = k 
    9    ppvalu = (ppvalu/fmmjdr)*h + coef(m,i)
         m = m - 1
         fmmjdr = fmmjdr - 1.
         if (fmmjdr .gt. 0.)            go to 9 
   99                                   return
      end function ppvalu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine interv ( xt, lxt, x, left, mflag )
!  from  * a practical guide to splines *  by C. de Boor    
!computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
!
!******  i n p u t  ******
!  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
!  lxt.....number of terms in the sequence  xt .
!  x.....the point whose location with respect to the sequence  xt  is
!        to be determined.
!
!******  o u t p u t  ******
!  left, mflag.....both integers, whose value is
!
!   1     -1      if               x .lt.  xt(1)
!   i      0      if   xt(i)  .le. x .lt. xt(i+1)
!   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
!   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
!
!        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
!        indicates that  x  lies outside the CLOSED interval
!        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
!        intervals is due to the decision to make all pp functions cont-
!        inuous from the right, but, by returning  mflag = 0  even if
!        x = xt(lxt), there is the option of having the computed pp function
!        continuous from the left at  xt(lxt) .
!
!******  m e t h o d  ******
!  The program is designed to be efficient in the common situation that
!  it is called repeatedly, with  x  taken from an increasing or decrea-
!  sing sequence. This will happen, e.g., when a pp function is to be
!  graphed. The first guess for  left  is therefore taken to be the val-
!  ue returned at the previous call and stored in the  l o c a l  varia-
!  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
!  essary since the present call may have nothing to do with the previ-
!  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!  ilo  and are done after just three comparisons.
!     Otherwise, we repeatedly double the difference  istep = ihi - ilo
!  while also moving  ilo  and  ihi  in the direction of  x , until
!                      xt(ilo) .le. x .lt. xt(ihi) ,
!  after which we use bisection to get, in addition, ilo+1 = ihi .
!  left = ilo  is then returned.
!
      implicit none
      integer left,lxt,mflag,   ihi,ilo,istep,middle
      real (kind=8) x,xt(lxt)
      data ilo /1/
      save ilo  
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
!
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
!
!              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
!              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
!
!           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
!     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
!**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
      go to 111
      end subroutine interv


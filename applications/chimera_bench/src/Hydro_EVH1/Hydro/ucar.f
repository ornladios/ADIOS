
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c-----------------------------------------------------------------------
c
c   ALF SUBROUTINES START HERE.
c
c-----------------------------------------------------------------------
c
c   Associated Legendre Function Computation Subroutines (v121305)
c
c   File: alf_sr_v121305
c
c-----------------------------------------------------------------------
c
c   Integrated Value Computations
c
c-----------------------------------------------------------------------
c
c
c
      subroutine iseries(u1sq,u2sq,n,pmm_1,pmm_2,spimm)
c-----------------------------------------------------------------------
c
c   Gerstl (1980, p.193): Identity to compute the required number of
c                         terms for a specified level of precision.
c
c   Note: For small integration regions, cancelling effects (loss of
c         significant figures) govern the maximum precision attainable,
c         NOT the number of terms.
c
c         For high degrees (Gerstl, 1980, p.193) must be determined
c         by the northern-most boundary, NOT an average of the north
c         and south.
c
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      logerror = dlog(1.d-16)
      usq = dmax1(u1sq,u2sq)
      n1 = n+1
      x0 = (1.d0+logerror)/dlog(usq)
      im = dint(x0) + 1
c-----------------------------------------------------------------------
c
c   Paul (1978, eq.25): Series expansion for sectorial integrals.
c
c-----------------------------------------------------------------------
      x = dble(2*im)
      xpm = dble(n) + x
      top = x - 3.d0
      den = x - 2.d0
      sum2 = 0.d0
      sum1 = 0.d0
      do i = 1, im-1
         xpmi = 1.d0/xpm
        ratio = top/den
        sum2 = (sum2 + xpmi)*ratio*u2sq
        sum1 = (sum1 + xpmi)*ratio*u1sq
        xpm = xpm - 2.d0
        top = top - 2.d0
        den = den - 2.d0
      enddo  !  i
      sum2 = (sum2 + 1.d0/xpm)*pmm_2*u2sq
      sum1 = (sum1 + 1.d0/xpm)*pmm_1*u1sq
      spimm = dabs(sum1 - sum2)
c-----------------------------------------------------------------------
      return
      end
c
c
c
      subroutine alfisct(nrst,nrfn,thet1,thet2,nmax,gs,pmm,pimm,zero,
     &                   nmax0)
c-----------------------------------------------------------------------
c
c   This is a Fortran 77 subroutine to compute scaled definite integrals
c   of sectorial Associated Legendre Functions (ALFs), as well as point
c   values of sectorial ALFs for the northern and southern boundaries of
c   each integration region.
c
c   Definite integration is with respect to the cosine of co-latitude.
c
c   For input, two arrays hold the northern spherical polar distance
c   (in radians) of the northern ('thet2') and southern ('thet1')
c   boundaries of the required integrals. The arrays 'thet2(nr)' and
c   'thet1(nr)' contain, respectively, the northern and southern
c   boundary arguments for the parallel band 'nr'. Inputs 'nrst' and
c   'nrfn' denote, respectively, the northern-most and southern-most
c   parallel bands for which ALF integrals and point values are
c   computed.
c
c   Other input arguments are: The maximum order 'nmax' up to which
c   ALFs will be computed; a global scale factor 'gs'; and 'nmax0',
c   which is used to define the first dimension of the 2-D arrays
c   'pmm' and 'pimm' (see below).
c
c   Integrals and point values of sectorial ALFs are computed from
c   order zero to 'nmax', for all parallel bands from 'nrst' to 'nrfn'.
c
c   'zero' is an internal 1-D array that is not intended to pass values
c   between this subroutine and the calling program. It is included in
c   the calling statement to allow the calling program to define its
c   dimension (see below).
c
c   Sectorial ALF integrals are returned in the 2-D array 'pimm',
c   in which the 'pimm(m+1,nr)' element holds the sectorial ALF
c   integral of order 'm' and parallel band 'nr'. Sectorial point
c   values are returned in the 2-D array 'pmm', in which the
c   'pmm(m+1,nr)'and 'pmm(m+1,nr+1)' elements hold the sectorial ALF
c   point values for order 'm', and for, respectively, the northern
c   and southern boundary of the parallel band 'nr'. Note that this
c   output assumes that the parallel bands adjoin each other, so that
c   the southern boundary of parallel band 'nr' is also the northern
c   boundary of parallel band 'nr+1'.
c
c   External: This subroutine calls the subroutine 'iseries'.
c
c   The following arrays should be dimensioned in the calling program
c   as follows:
c
c   - thet1(maxr)
c   - thet2(maxr)
c   - pmm(nmax0+1,maxr+1)
c   - pimm(nmax0+1,maxr)
c   - zero(maxr+1)
c
c   where:
c
c   - 'maxr' should be set in the calling program such that 'nrfn' is
c     never larger than 'maxr', and
c   - 'nmax0' should be set in the calling program such that 'nmax' is
c     never larger than 'nmax0'.
c
c   Note that 'nmax' should never be larger than 2700.
c
c   This subroutine was originally called 'sectoral6'.
c
c-----------------------------------------------------------------------
c
c   Fwd/Rev sectorial decision rule: Gerstl (1980, Eq. 13)
c   Sectorial recursion over increasing order: Paul (1978, Eq. 22a)
c   Sectorial recursion over decreasing order: Gleason (1985, Eq. 4.23)
c
c   Paul MK, (1978). Recurrence relations for integrals of associated
c      Legendre functions, Bulletin Geodesique, 52, pp177-190.
c   Gerstl M, (1980). On the recursive computation of the integrals of
c      the associated Legendre functions, manuscripta geodaetica, 5,
c      pp181-199.
c   Gleason, DM (1985). Partial sums of Legendre series via Clenshaw
c      summation, manuscripta geodaetica, 10, pp115-130.
c
c-----------------------------------------------------------------------
c
c     ORIGINAL PROGRAM:                           SIMON HOLMES, OCT 2004
c     MODIFIED TO ZERO HIGH-ORDER SECTORIALS:     SIMON HOLMES, JUN 2005
c     MODIFIED:                                   SIMON HOLMES, AUG 2005
c     MODIFIED:                                   SIMON HOLMES, NOV 2005
c
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      integer*4 nrst,nrfn,nmax,nmax0
      real*8 gs
      parameter(nm=2700,nm1=nm+1)
      real*8 thet1(*),thet2(*)
      real*8 rt(0:2*nm+5),rti(0:2*nm+5),dbl(0:2*nm+5)
      real*8 a(nm1),b(nm1),c(nm1),d(nm1)
      real*8 tmp1(nm),tmp2(nm)
      real*8 pmm(nmax0+1,*),pimm(nmax0+1,*)
      integer*4 zero(*)
c-----------------------------------------------------------------------
      if (nmax.gt.nmax0.or.nmax0.gt.2700) then
        write(6,6002)
 6002   format(///5x,
     &'***      Error in s/r alfisct: nmax > nmax0 or nmax0 > 2700 ***',
     &        //5x,
     &'       ***            Execution Stopped             ***',///)
        stop
      endif  !  nmax,nmax0
c-----------------------------------------------------------------------
       nrst1 = nrst+1
       nrfn1 = nrfn+1
       nmax1 = nmax+1
       xmax1 = nmax1*1.d0
      nmax2p = 2*nmax+5
        pi = 4.d0*datan(1.d0)
       pip = pi + 1.d-10
       pis = pi - 1.d-10
       pi2 = 0.5d0 * pi
      pi2p = pi2 + 1.d-10
      pi2s = pi2 - 1.d-10
c-----------------------------------------------------------------------
      do n = 1, nmax2p
         rt(n) = dsqrt(dble(n))
        rti(n) = 1.d0/rt(n)
        dbl(n) = dble(n)
      enddo  !  n
       rt(0) = 0.d0
      rti(0) = 0.d0
      dbl(0) = 0.d0
      do n = 1, nmax
        ns1 = n-1
         n1 = n+1
         n2 = 2*n
        a(n1) = 1.d0/dbl(n2+2)
        b(n1) = rt(n)*rt(n2-1)*rt(n2+1)*rti(ns1)
        c(n1) = rt(n2+1)*rti(n2)
        d(n1) = 2.d0*rt(n1)*rti(n+2)*rti(n2+3)*rti(n2+5)
      enddo  !  n
c-----------------------------------------------------------------------
      small = dlog10(1.d-200)
      gslog = dlog10(gs)
      do nr = nrst, nrfn
        th2 = thet2(nr)
        th1 = thet1(nr)
         u2 = dsin(th2)
         u1 = dsin(th1)
        if (th2.lt.-1.d-10.or.th1.gt.pip.or.th1.lt.th2) then
	  write(6,*)"th2,th1,pip:",th2,thet2(1),nr,thet2(nr)
          write(6,6000)
          stop
        endif  !  th1;th2
        if (th2.lt.pi2s.and.th1.gt.pi2p) then
          write(6,6001)
          stop
        endif  !  th1;th2
        u = dsin(th2)
        if (dabs(u).lt.1.d-10) u = dsin(th1)
        if (u.eq.1.d0) then
          mmax1 = nmax1
        else
          testm1 = (small - gslog)/dlog10(u) + 1.d0
          if (testm1.gt.xmax1) testm1 = xmax1
          mmax1 = nint(testm1)
        endif  !  u
        zero(nr) = mmax1
      enddo  !  nr
c
      if (nrfn.ne.-1) then
        u = dsin(thet1(nrfn))
        if (dabs(u).lt.1.d-10) u = dsin(thet1(nrfn-1))
        if (u.eq.1.d0) then
          mmax1 = nmax1
        else
          testm1 = (small - gslog)/dlog10(u) + 1.d0
          if (testm1.gt.xmax1) testm1 = xmax1
          mmax1 = nint(testm1)
        endif  !  u
        zero(nrfn+1) = mmax1
      endif  !  nrfn
c-----------------------------------------------------------------------
      do nr = nrst, nrfn
        do n1 = 1, nmax1
           pmm(n1,nr) = 0.d0
          pimm(n1,nr) = 0.d0
        enddo  !  n1
      enddo  !  nr
      if (nrfn.ne.-1) then
        do n1 = 1, nmax1
          pmm(n1,nrfn+1) = 0.d0
        enddo  !  n1
      endif  !  nrfn
c-----------------------------------------------------------------------
      if (nrst.ne.0) then
        th2 = thet2(nrst)
         u2 = dsin(th2)
        if (dabs(u2).lt.1.d-10) u2 = 0.d0
        pmm(1,nrst) = gs
        if (nmax.gt.0) pmm(2,nrst) = rt(3)*u2*gs
        mmax1 = zero(nrst)
        do n1 = 3, mmax1
          pmm(n1,nrst) = c(n1)*pmm(n1-1,nrst)*u2
        enddo  !  n1
      endif  !  nrst
      do nr1 = nrst1, nrfn1
        mmax1 = zero(nr1)
        th1 = thet1(nr1-1)
         u1 = dsin(th1)
        if (dabs(u1).lt.1.d-10) u1 = 0.d0
        pmm(1,nr1) = gs
        if (nmax.gt.0) pmm(2,nr1) = rt(3)*u1*gs
        do n1 = 3, mmax1
          pmm(n1,nr1) = c(n1)*pmm(n1-1,nr1)*u1
        enddo  !  n1
      enddo  !  nr1
c-----------------------------------------------------------------------
      do nr = nrst, nrfn
        th1 = thet1(nr)
        th2 = thet2(nr)
         u1 = dsin(th1)
         u2 = dsin(th2)
         t1 = dcos(th1)
         t2 = dcos(th2)
        if (dabs(t1).lt.1.d-10) then
          t1 =  0.d0
          u1 =  1.d0
        elseif (dabs(t2).lt.1.d-10) then
          t2 =  0.d0
          u2 =  1.d0
        elseif (dabs(u1).lt.1.d-10) then
          u1 =  0.d0
          t1 = -1.d0
        elseif (dabs(u2).lt.1.d-10) then
          u2 =  0.d0
          t2 =  1.d0
        endif
        u1sq = u1*u1
        u2sq = u2*u2
        mmax1 = zero(nr)
        if (u1.lt.u2.and.dabs(u1).gt.1.d-10) mmax1 = zero(nr+1)
        mmax  = mmax1-1
c
        pimm(1,nr) = gs*
     &               2.d0*dsin((th1-th2)*.5d0)*dsin(((th1+th2)*.5d0))
        usq = dmin1(u2sq,u1sq)
        if (usq.eq.0.d0) then
          xk = 2
        else
          xk = (mmax)/(mmax1*usq)
        endif  !  usq
        if (u1.eq.1.d0.or.u2.eq.1.d0) then
          xk = -1.d0
        endif  !  u1;u2
c
        if (xk.lt.1) then
          if (nmax.ge.1)
     &    pimm(2,nr) = rt(3)*.5d0*(t2*u2-th2-t1*u1+th1)*gs
          if (nmax.ge.2)
     &    pimm(3,nr) = (rt(3)*rt(5)*pimm(1,nr)+t2*pmm(3,nr)
     &                                   -t1*pmm(3,nr+1))*rti(3)*rti(3)
          do n1 = 4, mmax1
            pimm(n1,nr) = a(n1)*( b(n1)*pimm(n1-2,nr)
     &                           +2.d0*(t2*pmm(n1,nr)-t1*pmm(n1,nr+1)))
          enddo  !  n1
c
        else  !  xk
          pmm_2 = pmm(mmax1,nr  )
          pmm_1 = pmm(mmax1,nr+1)
          call iseries(u1sq,u2sq,mmax,pmm_1,pmm_2,spimm)
          pimm(mmax1,nr) = spimm
          if (mmax.ne.0) then
            pmm_2 = pmm(mmax,nr  )
            pmm_1 = pmm(mmax,nr+1)
          endif  !  mmax
          call iseries(u1sq,u2sq,mmax-1,pmm_1,pmm_2,spimm)
          if (nmax.ne.0) pimm(nmax,nr) = spimm
          do n1 = 2, mmax-1
            tmp1(n1) = d(n1)*dbl(n1+2)
            tmp2(n1) = d(n1)*(t1*pmm(n1+2,nr+1) - t2*pmm(n1+2,nr))
          enddo  !  n1
          do n1 = mmax-1, 2, -1
            pimm(n1,nr) = tmp1(n1)*pimm(n1+2,nr) + tmp2(n1)
          enddo  !  n1
        endif  !  xk
      enddo  !  nr
c-----------------------------------------------------------------------
 6000 format(///5x,
     &'***  Error in s/r alfisct: check that 0<=thet2<=thet1<=pi   ***',
     &        //5x,
     &'       ***            Execution Stopped             ***',///)
c-----------------------------------------------------------------------
 6001 format(///5x,
     &'***  Error in s/r alfisct:    integrals across the equator  ***',
     &        //5x,
     &'       ***            Execution Stopped             ***',///)
c-----------------------------------------------------------------------
      return
      end
c
c
c
      subroutine alfiord(nrst,nrfn,thet1,thet2,m,nmax,pmm,pimm,pint,
     &                   p1,p2,nmax0)
c-----------------------------------------------------------------------
c
c   This is a Fortran 77 subroutine to compute scaled, definite
c   integrals of non-sectorial Associated Legendre Functions (ALFs), as
c   well as point values of non-sectorial ALFs for the northern and
c   southern boundaries of each integration region, for a given order
c   'm' and for all degrees from 'm' to 'nmax'.
c
c   Definite integration is with respect to the cosine of co-latitude.
c
c   For input, two arrays hold the northern spherical polar distance
c   (in radians) of the northern ('thet2') and southern ('thet1')
c   boundaries of the required integrals. The arrays 'thet2(nr)' and
c   'thet1(nr)' contain, respectively, the northern and southern
c   boundary arguments for the parallel band 'nr'. Inputs 'nrst' and
c   'nrfn' denote, respectively, the northern-most and southern-most
c   parallel bands for which ALF integrals and point values are to be
c   computed.
c
c   Sectorial ALF integrals are input through the 2-D array 'pimm',
c   in which the 'pimm(m+1,nr)' element holds the sectorial ALF
c   integral of order 'm' and parallel band 'nr'. Sectorial ALF point
c   values are input through the 2-D array 'pmm', in which the
c   'pmm(m+1,nr)' and 'pmm(m+1,nr+1)' elements hold the sectorial ALF
c   point values for order 'm', and for, respectively, the northern and
c   southern boundary of the parallel band 'nr'. Note that this input
c   assumes that the parallel bands adjoin each other, so that the
c   southern boundary of parallel band 'nr' is also the northern
c   boundary of parallel band 'nr+1'.
c
c   Other input arguments are: The order 'm' for which the ALF
c   integrals and point values will be computed; the maximum degree
c   of the ALFs 'nmax'; and 'nmax0', which is used to define the first
c   dimension of the 2-D arrays 'pmm', 'pimm', 'p1' ,'p2' and 'pint'
c   (see below).
c
c   For a given order 'm', the subroutine returns all ALF integrals,
c   from degree 'm' to 'nmax', for all parallel bands from 'nrst' to
c   'nrfn'.
c
c   ALF integrals are returned in the 2-D array 'pint', in which the
c   'pint(n+1,nr)' element holds the integral of degree 'n' and
c   parallel band 'nr', for the specified order 'm'. ALF point values
c   are stored in the 2-D arrays 'p2' (northern boundary of 'nr') and
c   'p1' (southern boundary of 'nr'), in which the '(n+1,nr)' element
c   holds the point value of degree 'n' and parallel band 'nr' in both
c   cases (again for the specified order 'm').
c
c   Note that any global scaling of the computed ALFs enters through
c   the pre-computed (input) sectorial values in 'pmm' and 'pimm'.
c
c   The following arrays should be dimensioned in the calling program
c   as follows:
c
c   - thet1(maxr)
c   - thet2(maxr)
c   - pmm(nmax0+1,maxr+1)
c   - pimm(nmax0+1,maxr)
c   - pint(nmax0+2,maxr)
c   - p1(nmax0+2,maxr)
c   - p2(nmax0+2,maxr)
c
c   where:
c
c   - 'maxr' should be set in the calling program such that 'nrfn' is
c     never larger than 'maxr', and
c   - 'nmax0' should be set in the calling program such that 'nmax' is
c     never larger than 'nmax0'.
c
c   Note that 'nmax' should never be larger than 2700.
c
c   This subroutine was originally called 'pintm7'.
c
c-----------------------------------------------------------------------
c
c   Zonal/Tes'l recursion over increasing degree: Paul (1978, Eq. 20a)
c
c   Paul MK, (1978). Recurrence relations for integrals of associated
c      Legendre functions, Bulletin Geodesique, 52, pp177-190.
c
c-----------------------------------------------------------------------
c
c     ORIGINAL PROGRAM:                           SIMON HOLMES, SEP 2004
c     REMOVED DOUBLE COEFF ARRAYS:                SIMON HOLMES, JUL 2005
c     MODIFIED:                                   SIMON HOLMES, NOV 2005
c
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter(nm=2700,nm1=nm+1)
      real*8 thet1(*),thet2(*),pmm(nmax0+1,*),pimm(nmax0+1,*)
      real*8 rt(0:2*nm+5),rti(0:2*nm+5),dbl(0:2*nm+5)
      real*8 e(nm1),f(nm1),g(nm1),h(nm1)
      real*8 temp1(nm1),temp2(nm1),temp3(nm1)
      real*8 p1(nmax0+2,*),p2(nmax0+2,*),pint(nmax0+2,*)
c-----------------------------------------------------------------------
      save
      data ix,nmax_old/0,0/
       nmax1 = nmax+1
      nmax2p = 2*nmax+5
c-----------------------------------------------------------------------
      if (nmax.gt.nmax0.or.nmax0.gt.2700) then
        write(6,6002)
 6002   format(///5x,
     &'***      Error in s/r alfiord: nmax > nmax0 or nmax0 > 2700 ***',
     &        //5x,
     &'       ***            Execution Stopped             ***',///)
        stop
      endif  !  nmax,nmax0
c-----------------------------------------------------------------------
      if (m.gt.nmax) then
        write(6,6001)
 6001   format(///5x,'***  Error in s/r alfiord: m > Nmax  ***',
     $          //5x,'***        Execution Stopped         ***',///)
        stop
      endif  !  m
c-----------------------------------------------------------------------
      if (ix.eq.0.or.nmax.gt.nmax_old) then
        ix = 1
        nmax_old = nmax
        do n = 1, nmax2p
           rt(n) = dsqrt(dble(n))
          rti(n) = 1.d0/rt(n)
          dbl(n) = dble(n)
        enddo  !  n
         rt(0) = 0.d0
        rti(0) = 0.d0
        dbl(0) = 0.d0
      endif  !  ix
      m1 = m+1
      m2 = m+2
      m3 = m+3
      do n = m+2, nmax
        temp1(n+1) = rt(n*2+1)*rti(n+m)*rti(n-m)
        temp2(n+1) = rt(n*2+1)*rti(n+m)*rti(n-m)/dbl(n+1)
        temp3(n+1) = rt(n+m-1)*rt(n-m-1)*rti(n*2-3)
        e(n+1) =  temp1(n+1)*rt(n*2-1)
        f(n+1) =  temp1(n+1)*temp3(n+1)
        g(n+1) = -temp2(n+1)*rt(n*2-1)
        h(n+1) =  temp2(n+1)*temp3(n+1)*dbl(n-2)
      enddo  !  n
c-----------------------------------------------------------------------
      do nr = nrst, nrfn
        th1 = thet1(nr)
        th2 = thet2(nr)
         u1 = dsin(th1)
         u2 = dsin(th2)
         t1 = dcos(th1)
         t2 = dcos(th2)
        if (dabs(t1).lt.1.d-10) then
          t1 =  0.d0
          u1 =  1.d0
        elseif (dabs(t2).lt.1.d-10) then
          t2 =  0.d0
          u2 =  1.d0
        elseif (dabs(u1).lt.1.d-10) then
          u1 =  0.d0
          t1 = -1.d0
        elseif (dabs(u2).lt.1.d-10) then
          u2 =  0.d0
          t2 =  1.d0
        endif
        u1sq = u1*u1
        u2sq = u2*u2
        t1sq = t1*t1
        t2sq = t2*t2
          p1(m1,nr) =  pmm(m1,nr+1)
          p2(m1,nr) =  pmm(m1,nr  )
        pint(m1,nr) = pimm(m1,nr  )
        xrt = rt(m1*2+1)
        p1(m2,nr) = xrt*t1*p1(m1,nr)
        p2(m2,nr) = xrt*t2*p2(m1,nr)
        pint(m2,nr) = xrt/dbl(m2)*(u1sq*pmm(m1,nr+1)-u2sq*pmm(m1,nr))
        do n1 = m3, nmax1
          ex = e(n1)
          fx = f(n1)
            p1(n1,nr) = ex*t1*p1(n1-1,nr) - fx*p1(n1-2,nr)
            p2(n1,nr) = ex*t2*p2(n1-1,nr) - fx*p2(n1-2,nr)
          pint(n1,nr) = g(n1)*(u2sq*p2(n1-1,nr) - u1sq*p1(n1-1,nr))+
     &                  h(n1)*pint(n1-2,nr)
        enddo  !  n1
      enddo  !  nr
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c
c   ALF SUBROUTINES END HERE.
c
c-----------------------------------------------------------------------

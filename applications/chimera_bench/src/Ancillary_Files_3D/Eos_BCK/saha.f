c***********************************************************************
      subroutine saha
c***********************************************************************
      USE array_module
      USE eos_bck_module
      USE physcnst_module
      implicit double precision (a-h,o-z)
      logical test,first
      save
c                                                                      c
      parameter(xms=.7,xm0=2.,xm0ms=xm0-xms)
      parameter(third=1./3.,ba=-7.075,zero=0.0)
      parameter(a0=.067,c0=2.3798e-4)
      parameter(ncnvge = 10,epsf=4.e-7,epsg=4.e-8,epsh=4.e-7,epse=4.e-8)
      parameter (bsil = -8.448)
      parameter (one=1.,tenm3=1.e-3,dzamin=1.e-12)
c                                                                      c
      dimension w(2),x(2),y(2),z(2),e(5),f(5),gg(5),hh(5),a(4,4)
      dimension fvec(4),delta(4)
c                                                                      c
      data first /.true./
c***********************************************************************
      fpick(adum,bdum,xdum) =
     &     dim(sign(one,xdum),zero)*adum+dim(zero,sign(one,xdum))*bdum
c ******* x.lt.0:f=b *********   x.gt.0:f=a *****************x.eq.0:f=a*
      feps(j)   = e(j)**2/epse +f(j)**2/epsf +gg(j)**2/epsg
     * + hh(j)**2/epsh
      test(j)      = feps(j) .lt.1.0  .or.
     &             (abs(xhbck) .lt. 0.005 .and.
     &    4.* gg(j)**2 .lt. epsg .and. abs(hh(j))*ahbck/tbck .lt. 1.)
c************************************* if a guess or no nuclei *********
          if (first) then
      con        = 0.01
      tc         = 12.
      first      = .false.
          end if
      bad = .false.
      nfail = 0
      therm = c0 * tbck * sqrt(tbck) / dbck
      if( ((xabck .ne. zero) .or. (zabck .ne. zero)) .and.
     &   (theta .ne. zero) ) go to 2
1     call net
c********************************************* no nuclei needed ********
c      if( (xhbck .lt. con .or. nfail.ge.2) .and. (.not. bad)
c     & .and. (nfail.ne.0 .or. dbck.lt.1.e12*6.02e-16) ) return
      if( (xhbck .lt. con .or. nfail.ge.2) .and. (.not. bad) ) return
c******************************************choose guess and set up ****
              call guess
c*******************************************begin the iterations********
2      continue
c***********************************************************************
             do 10 iters = 1, ncnvge
              itsav = iters
              w(1) = theta
              x(1) = zabck
              y(1) = uhat
              z(1) = un
c*****************get the zero point************************************
c
                  call hvbub(e(1),f(1),gg(1),hh(1))
      if(fmbad) go to 11
      phi0 = phi
      u0 = upack
      ya0 = ya
      fm0 = fm
                  if (test(1)) go to 25
      dzap = dza
      duhatp = duhat
      dunp = dun
      dthetp = dtheta
      dza = 1.0/(8.0*ws*zabck)/xhbck
      duhat = sign(tenm3,gg(1))
      dun = sign(tenm3,hh(1)*0.5)
      dtheta = .001
                   if(iters.gt.1) then
      dza=sign(max(dzamin,min(abs(dza), .05*abs(dzap))),-dzap)
      duhat=sign(min(abs(duhat),.01*abs(duhatp)),-duhatp)
      dun=sign(min(abs(dun),.01*abs(dunp)),-dunp)
      dtheta=sign(min(abs(dtheta), .05*abs(dthetp)),-dthetp)
                   endif
      dza = sign(min(abs(dza),0.001*zabck),dza)
            x(2) = x(1) + dza
            zabck = x(2)
            call hvbub(e(2),f(2),gg(2),hh(2))
            if(fmbad) go to 11
            if(test(2)) go to 25
                  zabck = x(1)
            uhat = uhat + duhat
            y(2) = uhat
            call hvbub(e(3),f(3),gg(3),hh(3))
            if(fmbad) go to 11
            if(test(3)) go to 25
            uhat = y(1)
            un = un + dun
            z(2) = un
            call hvbub(e(4),f(4),gg(4),hh(4))
            if(fmbad) go to 11
            if(test(4)) go to 25
            un = z(1)
                 theta = theta + dtheta
                 w(2) = theta
                    call hvbub(e(5),f(5),gg(5),hh(5))
                    if(fmbad) go to 11
             if (test(5)) go to 25
c********************************************find derivs,fill matrix****
            a(1,1) = (f(2) - f(1))/(x(2) - x(1))
            a(1,2) = (f(3) - f(1))/(y(2) - y(1))
            a(1,3) = (f(4) - f(1))/(z(2) - z(1))
            a(1,4) = (f(5) - f(1))/(w(2) - w(1))
            a(2,1) = (gg(2) - gg(1))/(x(2) - x(1))
            a(2,2) = (gg(3) - gg(1))/(y(2) - y(1))
            a(2,3) = (gg(4) - gg(1))/(z(2) - z(1))
            a(2,4) = (gg(5) - gg(1))/(w(2) - w(1))
            a(3,1) = (hh(2) - hh(1))/(x(2) - x(1))
            a(3,2) = (hh(3) - hh(1))/(y(2) - y(1))
            a(3,3) = (hh(4) - hh(1))/(z(2) - z(1))
            a(3,4) = (hh(5) - hh(1))/(w(2) - w(1))
            a(4,1) = (e(2) - e(1))/(x(2) - x(1))
            a(4,2) = (e(3) - e(1))/(y(2) - y(1))
            a(4,3) = (e(4) - e(1))/(z(2) - z(1))
            a(4,4) = (e(5) - e(1))/(w(2) - w(1))
            fvec(1) = -f(1)
            fvec(2) = -gg(1)
            fvec(3) = -hh(1)
            fvec(4) = -e(1)
c********************************************invert the matrix**********
      call solve(4,a,fvec,delta,fmbad)
      if(fmbad) go to 11
      dza = delta(1)
      duhat = delta(2)
      dun = delta(3)
      dtheta = delta(4)
c***********************************limit step size*********************
      dza = max(yebck-.1-x(1),min(dza,.55-x(1)))
      dza = sign(min(abs(x(1)*0.25),abs(dza)),dza)
      dtheta = sign(min(theta*0.2,abs(dtheta)),dtheta)
c************************************* shuffle to keep upack .gt. 0 ****
      if( fm0 .gt. 1.0001 .and. iters .gt. 1) then
      unshuf = - fpick(tbck*(fm0-1.0)/(fm0+12.0*ya0),0.75*dzap ,-dzap)
      dunp = dun + unshuf
      dun = unshuf
      duhat = 1.0*tbck*(fm0-1.0)/(fm0+4.0*ya0)
      endif
c******************************** update the values ********************
      zabck = x(1) + dza
      uhat = y(1) + duhat
      un = z(1) + dun
      theta = w(1) + dtheta
c***********************************************************************
                if(iters.ge.5) then
      if((abs(dza/x(1))+abs(duhat/y(1))+abs(dun/z(1))).le.3.e-7)go to 20
      if(dbck .le. 0.1) go to 10
      if( (dbck .ge. d00*phi0) .or. (u0.gt.0.99)
     & .or. (iters.ge. ncnvge .and. bad) ) then
                        dtran = dbck
                        write(6,993)jshel,dtran,upack,d0,theta,xhbck
                        call eosnm
                        return
                 endif
                 endif
10                     continue
c******************************************************failure?*********
11    nfail = nfail + 1
      if(.not. bad) go to 1
      if(.not. fmbad) print *, jshel,dbck,tbck,yebck
      if(.not. fmbad) stop 'failure in Saha'
c******************************************************failure?*********
                        dtran = dbck
                        write(6,993)dtran,upack,d0,theta,xhbck
                        call eosnm
                        return
c******************************************************failure?*********
20    continue
25    continue
      xnbck = ynbck*u1pack
      xpbck = ypbck*u1pack
      xabck = 4.0*u1pack*ya
      pd  = dbck * (tsi*tbck)
      ed = 1.5*(tsi*tbck)*u1pack
      sd = tsi * u1pack * 2.5 - ( (1.-xhbck)*un 
     * -(yebck-zabck*xhbck)*uhat -xabck*ba)/tbck - xhbck * uhtra/tbck
      ph = xhbck * dbck * dbdlu
      sh = xhbck * af * tbck
     * * (2.*xmstar-2.*tbck/tc/(1.+tbck/tc)*(xmstar-xms) )
      eh = xhbck * b +xabck*ba + tbck*sh
      return
c***********************************************************************
993   format(' j=',i3,' set dtran=',f6.5,'  upack=',f7.3,'  d0=',f7.5,
     &' theta=',f7.4,' xhbck = ',1pe12.3)
999   format(' zero level for fe: yebck,d0,theta,wsize,egy',5f10.5)
1000  format(' i failed iters,zabck,uhat,un,theta ',/,2i3,4(1pe12.5))
                        end

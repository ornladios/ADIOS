      subroutine eos0
c***********************************************************************
      USE array_module
      USE eos_bck_module
      USE physcnst_module
      USE edit_module, ONLY : nlog
      implicit double precision (a-h,o-z)
      logical first
      save
c***********************************************************************
      parameter(third=1./3.,ba=-7.075,zero=0.0,dnp=1.2935)
      parameter(a0=.067,c0=2.3798e-4)
c                                                                      c
      data first/.true./
c***********************************************************************
      f3(z)        = sign( (abs(z))**third , z)
      fw1(x)       = 78.0 * x * x * (1.0-x) * f3(1.0-x)
      fbulk(x)     = wnm + ws * (1.0-2.0*x)*(1.0-2.0*x)
      fphi(x)      = 1.0 - 3.0 * (0.5-x)**2
      fxk(zabck)   = xk0*(1. - xkzafac*(1.- 2.*zabck)**2)
      fcomp(theta) = (fxk(zabck)/18.0) * (1.0 - theta)*(1.0 - theta)
c******************************************** get energy zeroes ********
              if(first)then
      d00 = 0.16
      yeFe = 26./56.
      zabck = yeFe
      theta0 = 1.0
      size0 = fw1(yeFe) /(fphi(yeFe))**third
      theta0 = 1 +(3.0/fxk(yeFe))*size0/theta0/theta0**third
      size0 = size0/theta0**third
      theta0 = 1 +(3.0/fxk(yeFe))*size0/theta0/theta0**third
      size0 = size0/theta0**third
      egy0 = fbulk(yeFe) + size0 + fcomp(theta0)
      write(nlog,999)yeFe,d00*theta0*fphi(yeFe),theta0,size0,egy0
      first = .false.
             endif
c***********************************************************************
                        return
c***********************************************************************
999   format(' zero level for fe: yebck,d0,theta,wsize,egy',5f10.5)
                        end

!///////////////////////////////////////////////////////////////////////

c***********************************************************************
      subroutine hvbub(tol1,tol2,tol3,tol4)
c***********************************************************************
      USE array_module
      USE eos_bck_module
      USE physcnst_module
      USE edit_module, ONLY : nlog
      implicit double precision (a-h,o-z)
      save
c                                                                      c
      parameter(zero=0.,one=1.,xms=.7,xm0=2.,xm0ms=xm0-xms,third=1./3.)
      parameter(bsil=-8.448,ba=-7.075)
      parameter(a0=.067,c0=2.3798e-4,piby2=1.570796327)
      common /tcrit/tc
c***********************************************************************
      f3(z) = sign( (abs(z))**third , z)
      frho(upack) = 1.0 -1.5*f3(upack) +0.5*upack + 1.e-8
      fg(upack,u1pack) = f3(frho(upack)) + f3(frho(u1pack))
      fw1(x) = 78.0 * x * x * (1.0-x) * f3(1.0-x)
      fwsize(zabck,upack,u1pack) = u1pack * fw1(zabck)
     *           * fg(upack,u1pack) / f3((phi*theta))
      fbulk(x)  = wnm + ws * (1.0-2.0*x)*(1.0-2.0*x)
      fdbldx(x) = - 4.0 * ws * (1.0-2.0*x)
      fdwszx(x) = wsize * 2.0 * (1.0-(5./3.)*x)/(x*(1.0-x))
      fphi(x) = 1.0 - 3.0 * (0.5-x)**2
      fder(upack) = ( f3(upack) - upack ) / f3(frho(upack))**2
      fdbdlu(upack) = -  fw1(zabck) / (phi*theta)**(third) *
     *              ( upack*fg(upack,u1pack) +(u1pack*fder(upack)
     *              -upack*fder(u1pack) )/6.0 )
      fah(zabck,upack) = 193.3333333 * (1.0-zabck)**2 / frho(upack)
      fuhtra(yh) = (tbck/ahbck) 
     *              *log (abs( yh/ (ahbck*sqrt(ahbck)*therm) ) )
      fxk(zabck) = xk0*(1. - xkzafac*(1. - 2.*zabck)**2)
      fcomp(theta) = (fxk(zabck)/18.0) * (1.0 - theta)*(1.0 - theta)
      fpcomp(theta,zabck) = -(fxk(zabck)/9.0)*theta*(1.0-theta)
      fbxthe(theta) = xk0*xkzafac*4./18.*(1. - 2.*zabck)
     *            *(1. - theta)**2
      faf(x) = sin(piby2*min(one,x)) / x**(2.*third)
      fcot(x) = cos(piby2*min(one,x)) / sin(piby2*min(one,x))
      fpth(x) = - af*xmstar*(-2.*third + piby2*x*fcot(x))
      fpick(adum,bdum,xdum) =
     &     dim(sign(one,xdum),zero)*adum+dim(zero,sign(one,xdum))*bdum
c****** x.lt.0:f=b *********   x.gt.0:f=a *****************x.eq.0:f=a*
c*********************************************************************
      tc              = 12.
      fmbad           = .false.
      if(abs(un/tbck) .ge. 500. .or. abs(uhat/tbck) .ge. 500.) then
      fmbad           = .true.
                      return
                      endif
      ynbck              = 2.0 * therm * exp((un/tbck))
      ypbck              = ynbck * exp(-(uhat/tbck))
      etaatr          = ( 4.0*(un-0.5*uhat-ba)/tbck )
      ya              = 8.0 * therm * exp(etaatr)
                      if( max(ynbck,ypbck,ya) .ge. 2.0) then
      fmbad           = .true.
                      return
                      endif
c*********************************************************************
      fm      = ynbck + ypbck + 4.0*ya
      phi     = fphi(zabck)
      dlphix  = 6.0*(0.5-zabck)/phi
      d0      = d00 * phi * theta
      upack       = (1.0 - fm) / (d0/dbck - fm)
      u1pack      = 1.0 - upack
      xhbck      = upack * d0/dbck
      ahbck      = fah(zabck,upack)
              if(ahbck .le. zero) then
      fmbad   = .true.
              return
              endif
              if ( u1pack .le. zero ) then
      xabck      = zero
      zabck      = zero
      dtran   = zero
      fmbad   = .true.
              return
              endif
      yh      = xhbck/(ahbck*u1pack)
      tsi     = ynbck + ypbck + ya + yh
c*********************************************************************
      wsize   = fwsize(zabck,upack,u1pack)
      uhtra   = fuhtra(yh)
      xmstar  = xms + xm0ms * wsize/size0 / (1.+tbck/tc)**2
      af      = a0 * faf(phi*theta)
      excit   = af*tbck*tbck
      tcorr   = 1.0 - excit*xm0ms/size0/(1.+tbck/tc)**2
      dbdlu   = fdbdlu(upack)*tcorr
     &        - 0.5*(uhtra+1.5*tbck/ahbck)
     *        * (f3(upack)-upack)/frho(upack)
      pcomp   = fpcomp(theta,zabck)
      dbdlth  = pcomp - third * wsize * tcorr 
     *        + fpth(theta*phi)*tbck*tbck
      dbdlph  = dbdlth - pcomp
      dbdx    = fdbldx(zabck) + fdwszx(zabck)*tcorr
     *        + (dbdlph-dbdlth)*dlphix
     &        + (uhtra+1.5*tbck/ahbck) *2.0/(1.0-zabck)
     &        + fbxthe(theta)
c***********************************************************************
      b       = fbulk(zabck) + wsize - excit*xmstar + fcomp(theta)
      tol1    = (dbdlu - dbdlth + tsi*tbck*dbck/d0)
     *        / (fxk(zabck)/9.*theta)
      tol2    = dbdx + uhat
      tol3    = yebck -  (2.0 * ya + ypbck)*u1pack - zabck * xhbck
      tol4    = -un +zabck*uhat + b + dbdlth + uhtra
c***********************************************************************
              return
c***********************************************************************
c***********************************************************************
      entry getzeros
c***********************************************************************
c******************************************** get energy zeroes ********
      yeFe = 26./56.
      zabck = yeFe
      theta0 = 1.0
      size0 = fw1(yeFe) /(fphi(yeFe))**third
      theta0 = 1 +(3.0/fxk(yeFe))*size0/theta0/theta0**third
      size0 = size0/theta0**third
      theta0 = 1 +(3.0/fxk(yeFe))*size0/theta0/theta0**third
      size0 = size0/theta0**third
      egy0 = fbulk(yeFe) + size0 + fcomp(theta0)
c***********************************************************************
      write(nlog,999)yeFe,d00*theta0*fphi(yeFe),theta0,size0,egy0
999   format(' zero level for fe: yebck,d0,theta,wsize,egy',5f10.5)
      return
c***********************************************************************
      entry guess
c***********************************************************************
      rxhmax = xhbck*fpick( yeFe/yebck,(1.-yeFe)/(1.-yebck), 
     *        yeFe-yebck)
      zabck = rxhmax*yebck + (1. - rxhmax)*yeFe
      phi = fphi(zabck)
      theta = 1.0+fw1(zabck)*(3.0/fxk(zabck))*(1.0-4.0*dbck/d00/phi)
      d0 = d00 * theta * phi
      excit = a0/(theta*phi)**(2.0*third) * tbck * tbck
      excit = a0 * tbck * tbck
      tcorr = 1.0 - excit*xm0ms/size0
      upack = dbck/d0 * yebck/zabck
      u1pack = 1.0 - upack
      wsize = fwsize(zabck,upack,u1pack)
      ahbck = fah(zabck,upack)
      uhathv = -fdbldx(zabck) - fdwszx(zabck) * tcorr
      unhv   = zabck * uhathv +(fbulk(zabck) + wsize - excit)
     *       + fdbdlu(upack)*tcorr +  fuhtra(1.0/ahbck/u1pack)
     *       + upack * tbck/ahbck/u1pack
c***********************************************************************
       un = unhv*rxhmax + (1.-rxhmax)*(un-.001)
       uhat = uhathv*rxhmax + (1. - rxhmax)*uhat
c***********************************************************************
                        return
c***********************************************************************
c************************************************ nuclear matter eos ***
      entry eosnm
c***********************************************************************
                        ahbck = zero
                        pd = zero
                        sd = zero
                        ed = zero
                        xabck = zero
                        xnbck = zero
                        xpbck = zero
                        xhbck = 1.0
                        upack = 1.0
                        zabck = yebck
                        d0 = dbck
                        xmstar = xms
              theta = (d0/d00)/fphi(yebck)
              af    = a0 * faf(d0/d00)
              excit = xmstar * af * tbck * tbck
      if (theta.gt.1.0) then
              ecomp = (fxk(zabck)/9.0/gamhv) *
     &            ( (theta**(gamhv-1.)-gamhv)/(gamhv-1.) + 1./theta)
              pcomp = (fxk(zabck)/9.0/gamhv)
     &                 * (theta**(gamhv-1.) - 1.0/theta)
      else
              ecomp = fcomp(theta)
              pcomp = fpcomp(theta,zabck)
      endif
              b = fbulk(yebck) - excit + ecomp
              ph = dbck * (pcomp + fpth(d0/d00) *tbck*tbck)
              sh = 2.0 * excit/tbck
              eh = b + tbck * sh
              dlphix = 6.0*(0.5-yebck)/fphi(yebck)
              uhat = -fdbldx(zabck) + pcomp*dlphix
     &             - ecomp/fxk(zabck)*xk0*xkzafac*4.*(1. - 2.*zabck)
              un = b + yebck * uhat + ph/dbck
                        return
c***********************************************************************
      end

!///////////////////////////////////////////////////////////////////////

c***********************************************************************
      subroutine net
c***********************************************************************
      USE array_module
      USE eos_bck_module
      USE physcnst_module
      implicit double precision (a-h,o-z)
      save
c                                                                      c
      parameter(nuc=6,zero=0.0,const=2.3798e-4,ncnvge=20,dnp=1.2935)
      parameter (expmax = 75.,three9 = 0.999, five9 = 0.99999)
c                                                                      c
      dimension a(nuc),z(nuc),zza(nuc),zan(nuc),a52(nuc),bn(nuc)
      dimension x(nuc),betem(nuc)
c                                                                      c
      data x/nuc*0./
      data a/1,1,4,56,28,16/,z/1,0,2,26,14,8/
      data zan/1,0,.5,.4643,.5,.5/,zza/1,0,1.,12.071,7,4/
      data bn/0,0,-7.075,-8.791,-8.448,-7.976/
      data betem/0,0,-7.075,-8.791,-8.448,-7.976/
      data a52/2,2,32.,23467.7,4148.54,1024/
c*** 1=p 2=n 3=2He4 4=26Fe56 5=14Si28 6=8O16****************************
c***********************************************************************
      fexp(f) = exp( min( expmax, max( -expmax , f) ) )
c***********************************************************************
      nflunk      = 0
      bad         = .false.
      epsp        = 5000.
                  if(uhat.ne. zero) go to 8
5                 continue
                  if (nflunk .le. 1) then
                  if(yebck.lt. 0.4643) then
      un          = tbck*log(0.5*(1.-yebck/.4643)/therm)
      uhat        = (un-bn(4)-tbck/a(4)*log(yebck/.4643/therm/a52(4)))
     *             /.4643
                  else
      xhbck       = (1.0-yebck)/(.5357)
      xpbck       = 1.-xhbck
      uhat        = (bn(4)+tbck/a(4)*log(xhbck/a52(4)/therm)
     &                    - tbck*log(xpbck*0.5/therm))/.5357
      un          = uhat + tbck*log(xpbck*0.5/therm)
                  endif
                  else if(nflunk.eq.2) then
      uhat        = tbck*log((1.0-yebck)/yebck)
      un          = tbck*log((0.5/therm)*(1.0-yebck))
                  else if (nflunk .eq. 3 .and. yebck .ge. .50) then
      xsil        = min((1-yebck)*2.,three9)
      xpbck       = 1.-xsil
      uhat        = (bn(5)+tbck/a(5)*log(xsil/a52(5)/therm)
     &               -tbck*log(xpbck*0.5/therm))/.5
      un          = uhat + tbck*log(xpbck*0.5/therm)
                  else if (nflunk .eq. 3) then
      xsil        = min(yebck*2.,three9)
      xpbck       = 1.-xsil
      uhat        = (bn(5)+tbck/a(5)*log(xsil/a52(5)/therm)
     &               -tbck*log(xpbck*0.5/therm))/.5
      un          = uhat + tbck*log(xpbck*0.5/therm)
                  else if(nflunk .eq. 4) then
      xhbck       = min(2.*yebck,five9)
      xpbck       = 1. - xhbck
      uhat        = (bn(4)+tbck/a(4)*log(xhbck/a52(4)/therm)
     &              - tbck*log(xpbck*0.5/therm))/.5357
      un          = uhat + tbck*log(xpbck*0.5/therm)
                  endif
8                 continue
c***********************************************************************
c                  if(exfac.eq.zero)then
c      write(6,*)' give exfac'
c      read(5,*)exfac
c                   endif
      exfac = zero
      excited = exfac*.067 * 2.0 * tbck * tbck
      do 90 n = 4,nuc
              betem(n) = bn(n) - excited
90            continue
c***********************************************************************
      do 100 itrat = 1,ncnvge
              ones=     zero
              zas =     zero
              zs =      zero
              as =      zero
              zzas =    zero
              bes  =    zero
              tsi =     zero
c***********************************************************************
              do 10 n = 1,nuc
                x(n) = fexp(log(a52(n)*therm)
     &                 + (a(n)*un-z(n)*uhat-a(n)*betem(n))/tbck )
                 ones=     ones + x(n)
                 zas =     zas + x(n)*zan(n)
                 zzas =    zzas + x(n)*zza(n)
                 zs =      zs + x(n)*z(n)
                 as =      as + x(n)*a(n)
                 bes =     bes + x(n)*bn(n)
                 tsi =     tsi + x(n)/a(n)
10            continue
c***********************************************************************
              tolone = log(ones)
              tolza = log(zas/yebck)
              eps = abs(tolone) + abs(tolza)
              if(eps.le. 2.e-5) go to 200
c***********************************************************************
              epssav = min(abs(eps),epsp)
              if(epssav .lt. epsp) then
                 epsp = epssav
                 xps = x(1)
                 xns = x(2)
                 xas = x(3)
                 xhs = x(4) + x(5) + x(6)
                 uns = un
                 uhats = uhat
              endif
c***********************************************************************
              if (tolone.ge.75. .and. tolza.ge.75.)go to 101
c***********************************************************************
              rdenom = ((zs/ones*zs/zas - as/ones*zzas/zas))
              if(rdenom .eq. zero) go to 101
              rdenom = tbck/rdenom
              dun = (tolone*zzas/zas - tolza/ones*zs)*rdenom
              duhat = (zs/zas*tolone - as/ones*tolza) * rdenom
              un = un + dun
              uhat = uhat + duhat
100   continue
c***********************************************************************
101   nflunk = nflunk + 1
c     write(6,1000)jshel,nflunk,tbck,dbck/6.022e-16,yebck
c     write(6,1001)un,uhat,dun,duhat,tolone,tolza
      if(nflunk.le.4) go to 5
            bad = .true.
1000  format(' net fail -jshel nflunk tbck dbck yebck'
     & ,2i4,3(1pe11.3))
1001  format(' un uhat dun duhat tolone tolza',6(1pe11.3))
            xpbck = xps
            xnbck = xns
            xabck = xas
            xhbck = xhs
            un = uns
            uhat = uhats
            go to 201
200   continue
      xpbck = x(1)
      xnbck = x(2)
      xabck = x(3)
      xhbck = x(4) + x(5) + x(6)
c***********************************************************************
201   continue
      upack  = zero
      theta = zero
      ahbck = zero
      zabck = zero
      xhsum = zero
           do 202 n = 4,nuc
      ahbck = ahbck + x(n)*a(n)
      zabck = zabck + x(n)*zan(n)
      xhsum = xhsum + x(n)
202        continue
      ahbck = ahbck/xhsum
      zabck = zabck/xhsum
      ph = zero
      sh = xhsum * excited * 2.0 / tbck
c Change made on 10/18/05 by SWB
      b = bes
!      b = bes - x(3)*bn(3)
      eh = bes  + xhsum*excited
      pd = dbck * tbck * tsi
      ed = 1.5 * tbck * tsi
      sd = 2.5 * tsi - (un - yebck*uhat - bes + xhsum*excited) /tbck
                        return
c***********************************************************************
                        end

!///////////////////////////////////////////////////////////////////////

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
                        write(6,993)dtran,upack,d0,theta,xhbck
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
993   format(' set dtran=',f6.5,'  upack=',f7.3,'  d0=',f7.5,
     &' theta=',f7.4,' xhbck = ',1pe12.3)
999   format(' zero level for fe: yebck,d0,theta,wsize,egy',5f10.5)
1000  format(' i failed iters,zabck,uhat,un,theta ',/,2i3,4(1pe12.5))
                        end

!///////////////////////////////////////////////////////////////////////

c***********************************************************************
c      comments of september 7,1984
c       solves the matrix equation a * x = b where
c     a is an n by n matrix, b and x n-component vectors
c
c
c***********************************************************************
c***********************************************************************
      subroutine solve(n,a,bb,x,fmbad)
c***********************************************************************
      implicit double precision (a-h,o-z)
      logical fmbad
      save
      parameter (nsize = 4)
      dimension a(nsize,nsize),bb(nsize),x(nsize)
      dimension xl(nsize,nsize),u(nsize,nsize),y(nsize)
      do 10 m=1,n
10    xl(m,m)=1.
        do 100 k=1,n
        do 40 j=k,n
        u(k,j)=a(k,j)
        if(k.eq.1)go to 40
        do 310 m=1,k-1
310      u(k,j)=u(k,j)-xl(k,m)*u(m,j)
40    continue
41      if(k.eq.n)go to 100
        do 50 i=k+1,n
        xl(i,k)=a(i,k)
        if(k.eq.1) go to 50
        do 45 m=1,k-1
45      xl(i,k)=xl(i,k)-xl(i,m)*u(m,k)
      if(u(k,k) .eq. 0.) then
      fmbad = .true.
      return
      endif
50      xl(i,k)=xl(i,k)/u(k,k)
100     continue
      do 200 i=1,n
      y(i)=bb(i)
      if(i.eq.1) go to 200
      do 150 k=1,i-1
150   y(i)=y(i)-xl(i,k)*y(k)
200   continue
      do 300 ii=1,n
      i=n+1-ii
      x(i)=y(i)
      if(i.eq.n)go to 300
      do 250 k=i+1,n
250   x(i)=x(i)-u(i,k)*x(k)
300   x(i)=x(i)/u(i,i)
      return
      end
c***********************************************************************

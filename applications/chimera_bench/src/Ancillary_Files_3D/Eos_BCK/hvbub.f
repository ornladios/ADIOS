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

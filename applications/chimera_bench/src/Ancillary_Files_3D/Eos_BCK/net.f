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

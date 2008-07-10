c-------------------------------------------------------------------   
      SUBROUTINE cc_difcs(e,ep,costh,T,s,dm,mun,mup,mue,gv,ga,fwqt)
c-------------------------------------------------------------------       
c-    THIS ROUTINE CALCULATES DIFF CROSS SECTION /VOLUME
c-        FOR  CHARGED CURRENT REACTIONS
c     input: e = incoming neutrino energy
c            ep = outgoing neutrino energy
c            costh = cosine of the angle between in and out neutrino 
c                    momentum vectors
c            T = temperature
c            s = mass of the target particle
c            mun = neutron chemical potential (including rest mass)
c            mup = proton chemical potential (  "   )
c            mue = electron chemical potential
c            dm = neutron - proton mass difference
c            gv = vector coupling
c            ga = axial coupling
c            fwqt = (1/V) dsigma/(d(costh) d(ep))
c    UNITS: if all energies are input in MeV the output fwqt has 
c           dimensions of MeV^-1 m^-1
c-------------------------------------------------------------------       
      Implicit none 
      
      real*8 e,ep,costh,fwqt,qo,qu2,q,T,s,dm,gv,ga
      real*8 pi,pi2,G2,forfac,mun,mup,mue,muh,efac,beta
      parameter (pi=3.141592654, pi2=9.8696044, G2= 6.79d-10)
c 
c                 G2 in MeV^-5 * m^-1
c

      real*8 eminus,l,u,z,arg1,arg2,pd,fd
      real*8 I0, I1, I2, uli1, uli2, uli3, bli1, bli2, bli3
      real*8 zi1, zi2, zi3, I2U, I2D, impl, impt, impa, impva
      real*8 r1, r2, r3, A, B, C

c------------------------------
c--      KINEMATICAL FACTORS          
c------------------------------
      qo = e - ep
      qu2 = 2.d0*E*EP*(costh-1.d0)
      q  = dsqrt(E*E + EP*EP - 2.d0*E*EP*costh)
      efac = (mue - ep)/T
      muh = mun - mup

      if(efac.gt.20.d0.or.qu2.ge.0.d0)then
         fwqt = 0.d0
         return
      else 
         fd = 1.d0/(1.d0+dexp(efac))
      end if
c------------------------------
c     -      RESPONSE FUNCTIONS
c------------------------------
      
      beta = 1.d0 - 2.d0*s*dm/qu2
      
      eminus = 0.5d0*(-qo*beta + q*dsqrt(beta**2 - 4.d0*(s*s/qu2)))
        
      if(eminus.lt.s) eminus = s
   
      l = eminus/T 
      u = mun/T
      z = (qo+muh)/T
      
      arg1 = l-u
      arg2 = l-u+z
      
      if(((arg1.gt.20.).and.(arg2.gt.20.)).or.(z.lt.-25.))then
         
         I0 = 0.d0
         I1 = 0.d0
         I2 = 0.d0
         
      else
         
         if (dabs(z).lt.1.d-8) then
            
            call polylog(arg1,uli1,uli2,uli3)
            
            zi1 = -1.d0/(1.d0+dexp(-arg1))
            zi2 = dlog(1.d0+dexp(arg1))
            zi3 = -uli2
            
            pd = 1.d0/(1.d0-0.5d0*z)
            
         else
            
            call polylog(arg1,uli1,uli2,uli3)
            call polylog(arg2,bli1,bli2,bli3)
            
            zi1 = (uli1 - bli1)/z
            zi2 = (uli2 - bli2)/z
            zi3 = (uli3 - bli3)/z
            pd = z/(1.d0-dexp(-z))
            
         endif
         
         I0 = T*(1.d0 + zi1)
         I1 = T*T*(u - 0.5d0*z + (zi2 + l*zi1))
         I2 = T*T*T*( u*(u-z) + (pi2+z*z)/3.d0 - 2.d0*zi3 
     &        + l*(2.d0*zi2 + l*zi1) )
         
      endif

      if(q.gt.1.d-3)then
         
         forfac = 1.d0
c------------------------
c---  long vector:
c------------------------
         impl = qu2*(I2 + qo*I1 + 0.25d0*qu2*I0 )/(2.d0*pi*q**3)
c------------------------
c---  tran vector:
c------------------------
         impt = 0.5d0*impl + (s*s+0.5d0*qu2)*I0/(4.d0*pi*q)
c------------------------
c---  axial correction
c------------------------
         impa = s*s*I0/(2.d0*pi*q)
         impva = qu2*(qo*I0+2.d0*I1)/(8.d0*pi*q**3)
         
      else
         
         forfac = dsqrt(0.5*(costh-1.d0)/(E*EP))
         
         impl =  forfac*(I2 + qo*I1 + 0.25d0*qu2*I0 )/(2.d0*pi)
         
         impt = 0.5d0*impl + forfac*(s*s+0.5d0*qu2)*I0/(4.d0*pi)
         
         impa = forfac * s*s*I0/(2.d0*pi)
         
         impva = forfac*(qo*I0+2.d0*I1)/(8.d0*pi)
         
      end if
      
      R1  = (gv*gv+ga*ga)*(impt+impl)
      R2  = gv*gv*impt+ga*ga*(impt-impa)
      R3  = 2.d0*gv*ga*impva
      
      A = (2.d0*E*EP+0.5d0*qu2)/(q*q)
      B = E + EP
      
      if(q.gt.1.d-3)then
         A = (2.d0*E*EP+0.5d0*qu2)/(q*q)
         B = E + EP           
         fwqt = (G2/pi2)*(costh-1.d0)*pd*(A*R1 + R2 + B*R3)
      else
         A = (2.d0*E*EP+0.5d0*qu2)/(q*q)
         B = E + EP
         fwqt = (G2/pi2)* pd *(A*R1 + R2 + B*R3)
      end if
      
      fwqt = fwqt*fd*ep**2
      
      return
      end

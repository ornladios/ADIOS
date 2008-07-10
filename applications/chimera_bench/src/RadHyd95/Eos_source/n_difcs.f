c-------------------------------------------------------------------   
      SUBROUTINE N_DIFCS(E,EP,costh,T,s,mun,cv,ca,fwqt)
c-------------------------------------------------------------------       
c-    THIS ROUTINE CALCULATES DIFF CROSS SECTION /VOLUME
c-        FOR  NEUTRAL CURRENT REACTIONS
c     input: E = incoming neutrino energy
c            EP = outgoing neutrino energy
c            costh = cosine of the angle between in and out neutrino 
c                    momentum vectors
c            T = temperature
c            s = mass of the target particle
c            mun = chemical potential (including rest mass of the target)
c            cv = vector coupling
c            ca = axial coupling
c            fwqt = (1/V) dsigma/(d(costh) d(EP))
c    UNITS: if all energies are input in MeV the output fwqt has 
c           dimensions of MeV^-1 m^-1
c-------------------------------------------------------------------       
      Implicit none 
      
      double precision E,EP,costh,fwqt,qo,qu2,q,T,s,mun,cv,ca
      double precision pi,pi2,G2
      parameter (pi=3.141592654, pi2=9.8696044, G2= 6.8945d-10)
c 
c                 G2 in MeV^-5 * m^-1
c

      double precision eminus,l,u,z,arg1,arg2,pd
      double precision I0, I1, I2, uli1, uli2, uli3, bli1, bli2, bli3
      double precision zi1, zi2, zi3, impl, impt, impa, impva
      double precision r1, r2, r3, A, B

c------------------------------
c--      KINEMATICAL FACTORS          
c------------------------------
      qo = E - EP
      qu2 = 2.d0*E*EP*(costh-1.d0)
      q  = dsqrt(E*E + EP*EP - 2.d0*E*EP*costh)

c------------------------------
c-      RESPONSE FUNCTIONS
c------------------------------

        eminus = 0.5d0*(-qo + q*dsqrt(1.d0 - 4.d0*(s*s/qu2)))
        if(eminus.lt.s) eminus = s
   
        l = eminus/T 
        u = mun/T
        z = qo/T

        arg1 = l-u
        arg2 = l-u+z

        if(((arg1.gt.25.).and.(arg2.gt.25.)).or.(z.lt.-25.))then
           I0 = 0.d0
           I1 = 0.d0
           I2 = 0.d0
        else

           if (dabs(z).lt.1.d-3) then

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
     &                + l*(2.d0*zi2 + l*zi1) )
              
        endif

c------------------------   
c---  long vector:
c------------------------  
        impl = qu2*(I2 + qo*I1 + 0.25d0*qu2*I0 )/(2.d0*pi*q**3)
c------------------------  
c---   tran vector:
c------------------------  
        impt = 0.5d0*impl + (s*s+0.5d0*qu2)*I0/(4.d0*pi*q)
c------------------------  
c---   axial correction
c------------------------  
        impa = s*s*I0/(2.d0*pi*q)
        impva = qu2*(qo*I0+2.d0*I1)/(8.d0*pi*q**3)
        
        R1  = (cv*cv+ca*ca)*(impt+impl)
        R2  = cv*cv*impt+ca*ca*(impt-impa)
        R3  = 2.d0*cv*ca*impva

        A = (2.d0*E*EP+0.5d0*qu2)/(q*q)
        B = E + EP

        fwqt = (G2/pi2)*(costh-1.d0)*pd*(A*R1 + R2 + B*R3)
     *         *EP*EP
         
      return
      end


c--------------------------------------------------   
      SUBROUTINE POLYLOG(x,pl1,pl2,pl3)
c--------------------------------------------------
      implicit none
      double precision x,pl1,pl2,pl3
      integer i,iload,nt,j
      parameter (nt=1000)
      double precision xa(0:nt),pl2a(0:nt),pl3a(0:nt)
      double precision pi2,arg,a,b
      parameter (pi2=9.8696044d0)

      save iload,xa,pl2a,pl3a
      data iload/0/

      if (iload.eq.0) then
        xa(0)=0.d0
        pl2a(0)=0.d0
        pl3a(0)=0.d0
c       open(unit=20,file='polylog.tab',status='unknown')
        open(unit=20,file='../../MGFLD/Input_Data/polylog.tab',
     &     status='old',position="rewind")
        do i=1,1000
          read(20,*) xa(i),pl2a(i),pl3a(i)
        enddo
        close(20)
        iload=1
      endif

      pl1 = x+dlog(1.d0+dexp(-x))
    
      if (x.ge.4.6) then
        pl2 = -(0.5d0*x*x + pi2/6.d0)/(1.d0+dexp(-1.5*x))
        pl3 = -x*(x*x + pi2)/6.d0/(1.d0+dexp(-1.69*x))
      else
        arg=dexp(x) 
        j = int(10.d0*arg)
        a = 10.d0*arg - dble(j)
        b = 1.d0 - a
        pl2 = pl2a(j)*b+pl2a(j+1)*a
        pl3 = pl3a(j)*b+pl3a(j+1)*a
      endif

      return
      end

c------------------------------------------------------------------


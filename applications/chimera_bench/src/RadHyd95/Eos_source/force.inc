C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       FORCE.INC
C    TYPE:         INCLUDE FILE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         2/29/92
C
C    PURPOSE:      This include file contains the statement function
C                  definitions for the nuclear force expressions.  To
C                  change the nucleon-nucleon interaction only this file
C                  need be changed.  NOTE:  We have assumed that the
C                  interaction is only density dependent, and that the
C                  INTERACTION HAS NO EXPLICIT TEMPERATURE DEPENDENCE!
C                  To introduce a temperature dependent interaction will
C                  require modification of the temperature derivatives
C                  in the routines NUCEOS & ALFEOS.
C
C
C    CALL LINE:    INCLUDE 'FORCE.INC/LIST'
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C                          Nucleon-Nucleon interaction parameters
      DOUBLE PRECISION AA, BB, CC, DD, DD3, SCRDD3
C
C                          Compression modulus, symmetry energy,
C                          binding energy, & saturation density
      DOUBLE PRECISION COMP, SYMM, BIND_E, NSUBS
C
C                          Surface symmetry energy & surface tension
      DOUBLE PRECISION SYM_S, SIG_S
C
C                          Numerical coefficient (called alpha in
C                          the Nucl. Phys. A, vol. 535, pg. 331 paper)
      DOUBLE PRECISION SKYRMC
C
C
C                   This common block contains the variables that
C                   specify the particular interaction used
      COMMON /SKYRME/ AA, BB, CC, DD, DD3, SCRDD3,
     1                COMP, SYMM, SKYRMC, BIND_E,
     2                NSUBS, SYM_S, SIG_S
C
      DOUBLE PRECISION PV_E, DPVEDN, DPVEDX
      DOUBLE PRECISION PV_PR, DPVRDP, DPVRDN
      DOUBLE PRECISION PVP, DPVPDP, DPVPDN, DVP_DX, DVP_DI
      DOUBLE PRECISION PVN, DPVNDP, DPVNDN, DVN_DX, DVN_DI
      DOUBLE PRECISION DENOM, VETERM, VTERM1, VTERM2
      DOUBLE PRECISION PROT_D, NUT_D, DENSIT, NNN, XXX
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------
C                 These terms enter in many of the equations
C-----------------------------------------------------------------------
C
      DENOM(DENSIT)=1.0+DD3*(DENSIT**(DD-1.0D+00))
C
      VETERM(DENSIT)= (CC*(1.0+DD)*(DENSIT**DD)+
     1     2.0D+00*CC*DD3*(DENSIT**(2.0D+00*DD-1.0D+00)))/
     2     (DENOM(DENSIT)**2)
C
C-----------------------------------------------------------------------
C                        Internal energy stuff
C-----------------------------------------------------------------------
C
C                 Internal energy (IE) due to interaction
      PV_E(PROT_D,NUT_D)=AA*((PROT_D+NUT_D)**2)+
     1     4.0D+00*BB*PROT_D*NUT_D+( CC*((PROT_D+NUT_D)**(1.0D+00+DD))/
     2     DENOM(PROT_D+NUT_D) )+PROT_D*DELTAM
C
C                 Derivative of IE w.r.t. density at fixed X
      DPVEDN(NNN,XXX)=2.0D+00*NNN*
     1    (AA+4.0D+00*BB*XXX*(1.0D+00-XXX))+VETERM(NNN)+XXX*DELTAM
C
C                 Derivative of IE w.r.t. X at fixed density
      DPVEDX(NNN,XXX)=4.0D+00*BB*(1.0D+00-2.0D+00*XXX)*(NNN**2)+
     1    NNN*DELTAM
C
C-----------------------------------------------------------------------
C                        Interaction potential stuff
C-----------------------------------------------------------------------
C
C                 Proton & neutron interaction potentials
      PVP(PROT_D,NUT_D)=2.0D+00*AA*(PROT_D+NUT_D)+
     1     4.0D+00*BB*NUT_D+VETERM(PROT_D+NUT_D)+DELTAM
C
      PVN(PROT_D,NUT_D)=2.0D+00*AA*(PROT_D+NUT_D)+
     1     4.0D+00*BB*PROT_D+VETERM(PROT_D+NUT_D)
C
C                 These terms enter in all of the equations
      VTERM1(DENSIT)=CC*(
     1    DD*(1.0D+00+DD)*(DENSIT**(DD-1.0D+00))+
     2    2.0D+00*DD3*(2.0D+00*DD-1.0D0)*(DENSIT**(2.0D0*DD-2.0D0))
     3    )/(DENOM(DENSIT)**2)
      VTERM2(DENSIT)=-2.0D+00*CC*DD3*(
     1    (DD**2-1.0D+00)*(DENSIT**(2.0D+00*DD-2.0D+00))+
     2    2.0D+00*DD3*(DD-1.0D+00)*(DENSIT**(3.0D+00*DD-3.0D+00))
     3    )/(DENOM(DENSIT)**3)
C
C                 Derivative of proton potential w.r.t.
C                 proton density at fixed neutron density
      DPVPDP(PROT_D,NUT_D)=2.0D+00*AA+
     1    VTERM1(PROT_D+NUT_D)+VTERM2(PROT_D+NUT_D)
C
C                 Derivative of proton potential w.r.t.
C                 neutron density at fixed proton density
      DPVPDN(PROT_D,NUT_D)=2.0D+00*AA+4.0D+00*BB+
     1    VTERM1(PROT_D+NUT_D)+VTERM2(PROT_D+NUT_D)
C
C                 Derivative of proton potential w.r.t.
C                 X at fixed density
cc      DVP_DX(NNN,XXX)=-4.0D+00*BB*NNN
C
C                 Derivative of proton potential w.r.t.
C                 density at fixed X
cc      DVP_DI(NNN,XXX)=2.0D+00*AA+4.0D+00*BB*(1.0D+00-XXX)+
cc     1    VTERM1(NNN)+VTERM2(NNN)
C
C                 Derivative of neutron potential w.r.t.
C                 neutron density at fixed proton density
      DPVNDN(PROT_D,NUT_D)=2.0D+00*AA+
     1    VTERM1(PROT_D+NUT_D)+VTERM2(PROT_D+NUT_D)
C
C                 Derivative of neutron potential w.r.t.
C                 proton density at fixed neutron density
      DPVNDP(PROT_D,NUT_D)=2.0D+00*AA+4.0D+00*BB+
     1    VTERM1(PROT_D+NUT_D)+VTERM2(PROT_D+NUT_D)
C
C                 Derivative of neutron potential w.r.t.
C                 X at fixed density
cc      DVN_DX(NNN,XXX)=4.0D+00*BB*NNN
C
C                 Derivative of neutron potential w.r.t.
C                 density at fixed X
cc      DVN_DI(NNN,XXX)=2.0D+00*AA+4.0D+00*BB*XXX+
cc     1    VTERM1(NNN)+VTERM2(NNN)
C
C-----------------------------------------------------------------------
C                        Pressure potential stuff
C-----------------------------------------------------------------------
C
C                 Interaction contribution to the pressure (IP)
      PV_PR(PROT_D,NUT_D)=PROT_D*PVP(PROT_D,NUT_D)+
     1    NUT_D*PVN(PROT_D,NUT_D)-PV_E(PROT_D,NUT_D)
C
C                 Derivative of IP w.r.t. proton density at
C                 fixed neutron density
      DPVRDP(PROT_D,NUT_D)=PROT_D*
     1    DPVPDP(PROT_D,NUT_D)+NUT_D*DPVNDP(PROT_D,NUT_D)
C
C                 Derivative of IP w.r.t. neutron density at
C                 fixed proton density
      DPVRDN(PROT_D,NUT_D)=PROT_D*
     1    DPVPDN(PROT_D,NUT_D)+NUT_D*DPVNDN(PROT_D,NUT_D)
C
C

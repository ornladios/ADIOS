C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EL_REL.FOR
C
C***********************************************************************
C
C    MODULE:       EL_REL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         2/12/91
C
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: swesty@nuclear.physics.sunysb.edu or
C                            fswesty@sbast3.ess.sunysb.edu
C
C    PURPOSE:      The relativistic electron and photon EOS
C
C
C    CALL LINE:    CALL EL_REL(T,YE,BRYDNS)
C
C    INPUTS:       T = TEMPERATURE
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:      NONE
C
C
C
C    INCLUDE FILES:  EL_REL.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EL_REL(T,YE,BRYDNS)
C
      USE el_eos_module
C
      IMPLICIT NONE
C
      DOUBLE PRECISION T, YE, BRYDNS
C
C
C
C                           Plancks constant & speed of light
      DOUBLE PRECISION HBAR, C
      PARAMETER (HBAR=6.58217317D-22,C=2.997924581D23)
C
C                           Pi and 1/3
      DOUBLE PRECISION PI, PI2, OVR3, MOVR3, OVR23
      PARAMETER(PI=3.1415927, PI2=9.8696044)
      PARAMETER(OVR3=0.33333333, MOVR3=-0.33333333, OVR23=0.66666667)
C
C                           2nd Fermi integral
      DOUBLE PRECISION F_2
C
C                           Positron degeneracy parameter
      DOUBLE PRECISION ELPETA
C
C
C                    Leptons
C
C                    Electron number density
      NSUBE = BRYDNS*YE
C
C                    Coefficants for chemical potential
C                    and thermodynamics quantities
      QSUBE = 1.0/( 3.0*(PI**2)*((HBAR*C)**3) )
C
      ACOEF = 0.5*NSUBE/QSUBE
C
      BCOEF = (ACOEF**2+((PI**6)*T**6)/27.0)**0.5
C
      DBDT = (PI**6)*(T**5)/(9.0*BCOEF)
C
      CCOEF = (ACOEF+BCOEF)**OVR3
C
C
C                    Electron chemical potential
      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF
C
C                    Positron degeneracy parameter
      ELPETA = -MUSUBE/T
C
C                    Positron number density
      NEPLUS =  3.0*QSUBE*(T**3)*F_2(ELPETA)
C
C
C                    Electron pressure for rel. case
      EPRESS = 0.25*QSUBE*(MUSUBE**4+2.0*(PI*T*MUSUBE)**2+
     1 7.0*((PI*T)**4)/15.0)
C
C
C                    Electron internal energy per baryon
      EU = 0.75*QSUBE*(MUSUBE**4+2.0*(PI*MUSUBE*T)**2+
     1 7.0*((PI*T)**4)/15.0)/BRYDNS
C
C
C                    Electron free energy per baryon
      FSUBE = ((MUSUBE*NSUBE)-EPRESS)/BRYDNS
C
C                    Electron entropy per baryon
      ES = QSUBE*(((PI*MUSUBE)**2)*T+7.0*(PI**4)*(T**3)/
     1 15.0)/BRYDNS
C
C                    Photons
C
C                    Photon pressure
      PPRESS = (PI**2)*(T**4)/(45.0*((HBAR*C)**3))
C                    Photon entropy per baryon
      PS = 4.0*PPRESS/(T*BRYDNS)
C
C                    Photon internal energy per baryon
      PU = 3.0*PPRESS/BRYDNS
C
C                    Photon free energy per baryon
      PF = PU-T*PS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C                    Derivatives of chem. potential w.r.t. T,
C                    BRYDNS, YE
C
      DEMUDT = DBDT/(3.0*CCOEF**2)-OVR23*(PI**2)*T/CCOEF+
     1         DBDT*((PI*T)**2)/(9.0*CCOEF**4)
C
      DEMUDN = (YE*PI2*(HBAR*C)**3)/(MUSUBE**2+OVR3*PI2*T**2)
C
      DEMUDY = BRYDNS*DEMUDN/YE
C
C
C                    Derivatives of pressure w.r.t. BRYDNS,YE,T
C
      DEPDN = BRYDNS*YE*DEMUDN
C
      DEPDY = BRYDNS*DEPDN/YE
C
      DEPDT = BRYDNS*(ES+YE*DEMUDT)
C
C
C                    Derivatives of entropy w.r.t. T,BRYDNS,YE
C
      DESDT = ES/T+OVR23*(7.0*PI2*(T**2)/15.0+MUSUBE*T*DEMUDT)/
     1        (BRYDNS*(HBAR*C)**3)
C
      DESDN = -1.0*DEPDT/(BRYDNS**2)
C
      DESDY = 2.0*T*QSUBE*PI2*MUSUBE*DEMUDY/BRYDNS
C
C
C                    Derivatives of internal energy w.r.t.
C                    T,BRYDNS,YE
      DEUDT = T*DESDT
C
      DEUDN = (YE*(MUSUBE-T*DEMUDT)-EU)/BRYDNS
C
      DEUDY = 3.0*QSUBE*((MUSUBE**3)+PI2*(T**2)*MUSUBE)*
     1        DEMUDY/BRYDNS
C
C
C                               Photons
C
C                    Derivatives of photon pressure
      DPPDN = 0.0
      DPPDT = BRYDNS*PS
      DPPDY = 0.0
C
C                    Derivatives of photon entropy
      DPSDN = -PS/BRYDNS
      DPSDT = 3.0*PS/T
      DPSDY = 0.0
C
C                    Derivatives of internal energy
      DPUDN = -0.75*T*PS/BRYDNS
      DPUDT = 3.0*PS
      DPUDY = 0.0
C
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
 999  RETURN
C
C
      END

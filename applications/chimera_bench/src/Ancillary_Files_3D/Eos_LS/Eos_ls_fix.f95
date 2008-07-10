Calfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalf
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         ALFEOS.FOR
C
C***********************************************************************
C
C    MODULE:       ALFEOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         8/30/90 Modified from model 4-A
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU
C                            FSWESTY@SBAST3.SUNYSB.EDU
C
C
C    CALL LINE:    CALL ALFEOS(INPVAR,YE,BRYDNS)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:
C
C
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE ALFEOS(INPVARP,YEP,BRYDNSP,P_PREVP,SSFLAGP)
C
C
      USE el_eos_module
      USE eos_m4c_module
C
      IMPLICIT NONE
C
      INTEGER SSFLAGP
      DOUBLE PRECISION INPVARP(4),YEP,BRYDNSP,P_PREVP
C
C                       "ZERO" flag
      INTEGER ZFLAG
C
C                       "Negative" flag
      INTEGER NFLAG
C
C                       Function type declarations
C
      DOUBLE PRECISION F_1_2, F_3_2, FINV12, FHALFO, FHALF
C
C                       Include the nucleon-nucleon interaction
C                       statement function definitions
      INCLUDE 'force.inc'
C
C
      INPVAR(1)      = INPVARP(1)
      INPVAR(2)      = INPVARP(2)
      INPVAR(3)      = INPVARP(3)
      INPVAR(4)      = INPVARP(4)
      YE             = YEP
      BRYDNS         = BRYDNSP
      P_PREV         = P_PREVP
      SSFLAG         = SSFLAGP
C
C
C
C                         Unset the zero flag
      ZFLAG = 0
C
C                         Ratio of baryon density to saturation density
      Y = BRYDNS/NSUBS
C
C
C                         Set T equal to the input variable (the entropy
C                         and internal energy options are not implemented
C                         in this version)
      T = INPVAR(1)
C
C
C                         Calc the quantum concentration of nucleons
      NQ = 2.36D-4*T**1.5
C
C                         Calc the Fermi integral coefficent
      UQ = 20.721
C
      MQ = (T/UQ)**1.5
C
      KQ = ((T/UQ)**2.5)/(2.0*PI**2)
C
      LQ = UQ*(MQ**OVR53)/(3.0*(PI**2))
C
      ETAMAX = 0.95*FINV12(2.0*(PI**2)*BRYDNS/MQ)
C
C
C
C                              Set the proton density to its old value
      NPOUT = P_PREV
C
      IF(BRYDNS.GT.(0.98*2.0/(YE*V_ALFA))) THEN
        NPOUT = YE*BRYDNS
        NNOUT = (1.0-YE)*BRYDNS
        NOUT = BRYDNS
C
C
c20        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)
        VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
C
c20        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+
c20     1    CC*(1.0+DD)*NOUT**DD+DELTAM)
        VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
        ZNO = 2.0*(PI**2)*NNOUT/MQ
C
        ZPO = 2.0*(PI**2)*NPOUT/MQ
C
        ETA_NO = FINV12(ZNO)
C
        ETA_PO = FINV12(ZPO)
C
        F32_NO = F_3_2(ETA_NO)
        F32_PO = F_3_2(ETA_PO)
C
        TAU_NO = KQ*F32_NO
        TAU_PO = KQ*F32_PO
C
C
c20        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(
c20     1    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
C
        MUN_O = T*ETA_NO+VNOUT
        MUN = MUN_O
C
        MUP_O = T*ETA_PO+VPOUT
        MUPROT = MUP_O
C
C                              Calculate diff. of chem. potentials
        MUHAT = MUN-MUPROT
C
C
C                              Calculate the alpha particle
C                              chemical potential
        MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA
C
        ALFDNS = 0.0
C
        EXALFA = 1.0-ALFDNS*V_ALFA
C
      ELSE
C
C                              Calculate the neutron density
        NNOUT = 2.0*BRYDNS*(1.0-2.0*YE)/(2.0-BRYDNS*YE*V_ALFA)+
     1            NPOUT*(2.0-(1.0-YE)*BRYDNS*V_ALFA)/
     2            (2.0-BRYDNS*YE*V_ALFA)
C
C                              Calculate density of outside nucleons
        NOUT = NPOUT+NNOUT
C
c20        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)
        VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
C
c20        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+
c20     1    CC*(1.0+DD)*NOUT**DD+DELTAM)
        VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
        ZNO = 2.0*(PI**2)*NNOUT/MQ
C
        ZPO = 2.0*(PI**2)*NPOUT/MQ
C
        ETA_NO = FINV12(ZNO)
C
        ETA_PO = FINV12(ZPO)
C
        F32_NO = F_3_2(ETA_NO)
        F32_PO = F_3_2(ETA_PO)
C
        TAU_NO = KQ*F32_NO
        TAU_PO = KQ*F32_PO
C
C
c20        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(
c20     1    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
C
        MUN_O = T*ETA_NO+VNOUT
        MUN = MUN_O
C
        MUP_O = T*ETA_PO+VPOUT
        MUPROT = MUP_O
C
C                              Calculate diff. of chem. potentials
        MUHAT = MUN-MUPROT
C
C
C                              Calculate the alpha particle
C                              chemical potential
        MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA
C
C                              Calculate density of alpha particles
C
        IF(ABS(MUALFA/T).LT.30.0) THEN
          ALFDNS = 8.0*NQ*DEXP(MUALFA/T)
        ELSEIF((MUALFA/T).LT.-30.0) THEN
          ALFDNS = 0.0
        ELSE
          ALFDNS = 8.0*NQ*DEXP(3.0D1)
        ENDIF
C
C
        EXALFA = 1.0-ALFDNS*V_ALFA
C
C                              Calculate "non-zeroness" of baryon
C                              conservation equation and save the
C                              value to be used in the finite
C                              difference approximation of DGDPRT
        GOLD = BRYDNS-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS
        PRTOLD = NPOUT
C
C                              Take a small step to get derivative
        NPOUT = NPOUT+0.001*BRYDNS
C
        DO 11 I=1,30,1
C
C                              Unset the negative flag
          NFLAG = 0
C
 14       CONTINUE
C                              Calculate the neutron density
          NNOUT = 2.0*BRYDNS*(1.0-2.0*YE)/(2.0-BRYDNS*YE*V_ALFA)+
     1            NPOUT*(2.0-(1.0-YE)*BRYDNS*V_ALFA)/
     2            (2.0-BRYDNS*YE*V_ALFA)
C
          IF((NNOUT.LT.0.0).AND.(I.EQ.1)) THEN
            NPOUT = PRTOLD-0.5*DPRT
          ELSEIF((NNOUT.LT.0.0).AND.(I.EQ.1).AND.(NFLAG.NE.1)) THEN
            NPOUT = 0.99*P_PREV
            NFLAG = 1
          ELSEIF((NNOUT.LT.0.0).AND.(I.EQ.1).AND.(NFLAG.EQ.1)) THEN
            SSFLAG = 0
            GOTO 999
          ENDIF
C                              Calculate density of outside nucleons
          NOUT = NPOUT+NNOUT
C
c20          VNOUT = EIFLAG*
c20     1      (2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)
          VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
C
c20          VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+
c20     1      CC*(1.0+DD)*NOUT**DD+DELTAM)
          VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
          ZNO = 2.0*(PI**2)*NNOUT/MQ
C
          ZPO = 2.0*(PI**2)*NPOUT/MQ
C
          ETA_NO = FINV12(ZNO)
C
          ETA_PO = FINV12(ZPO)
C
          F32_NO = F_3_2(ETA_NO)
C
          F32_PO = F_3_2(ETA_PO)
C
          TAU_NO = KQ*F32_NO
          TAU_PO = KQ*F32_PO
C
c20          BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*
c20     1      (AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )
          BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
          MUN_O = T*ETA_NO+VNOUT
          MUN = MUN_O
C
          MUP_O = T*ETA_PO+VPOUT
          MUPROT = MUP_O
C
C                              Calc difference of potentials
          MUHAT = MUN-MUPROT
C
C                              Calc alpha particle chemical potentials
          MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA
C
C                              Calc alpha particle density
C
          IF(ABS(MUALFA/T).LT.30.0) THEN
            ALFDNS = 8.0*NQ*DEXP(MUALFA/T)
          ELSEIF((MUALFA/T).LT.-30.0) THEN
            ALFDNS = 0.0
          ELSE
            ALFDNS = 8.0*NQ*DEXP(3.0D1)
          ENDIF
C
C
          EXALFA = 1.0-ALFDNS*V_ALFA
C
C                              Calc "non-zeroness" of baryon cons. eqn.
          G = BRYDNS-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS
C
C                              Calculate derivative of baryon conservation
C                              equation w.r.t. proton density by finite
C                              diference approximation
          DGDPRT = (G-GOLD)/(NPOUT-PRTOLD)
C
C                              If rate of change is near zero
C                              and zero flag is not set
          IF((ABS(DGDPRT).LT.1.0D-25).AND.(ZFLAG.EQ.0)) THEN
C                              Tweak the step size
            NPOUT = PRTOLD-0.5*DPRT
C                              and set the zero flag
            ZFLAG = 1
C                              and go back and try again
            GOTO 14
C                              If failure occurs again
          ELSEIF((ABS(DGDPRT).LT.1.0D-25).AND.(ZFLAG.EQ.1)) THEN
C                              declare an EOS failure
            SSFLAG = 0
C                              and return
            GOTO 999
          ENDIF
C
C                              Calculate new Newton-Raphson step
          DPRT = G/DGDPRT
C
C                              Save old value of proton density & G
          PRTOLD = NPOUT
          GOLD = G
C
C
 13       CONTINUE
C
C                              Potential "new" value of proton density
          PRTNEW = NPOUT-DPRT
C
C                              If new proton density is less than the
C                              baryon density and greater than zero
C                              then update the proton density
          IF(PRTNEW*(BRYDNS-PRTNEW).GT.0.0) THEN
            NPOUT = NPOUT-DPRT
C                              Else cut the step size in half and try again
          ELSE
            DPRT = DPRT*0.5
            GOTO 13
          ENDIF
C
C                              If step size is small enough break out of
C                              the DO 11 loop, otherwise continue
          IF(ABS(DPRT/NPOUT).LT.10E-11) GOTO 12
 11     CONTINUE
C
c      write(*,*) 'A failed to converge; switching to F' ! take out later
        SSFLAG = 0
        GOTO 999
C
C
 12     CONTINUE
C
      ENDIF
C                              Set the success flag
      SSFLAG = 1
C
C
C                              Calc outside nucleon density
      NOUT = NNOUT+NPOUT
C
C                              Calc outside nucleon fraction
      XOUT = NPOUT/NOUT
C
C                              Calculate particle fractions
      XALFA = 4.0*ALFDNS/BRYDNS
      XPROT = EXALFA*NPOUT/BRYDNS
      XNUT = EXALFA*NNOUT/BRYDNS
      XH = 0.0
C
C                              Baryons
C
      F32_NO = F_3_2(ETA_NO)
C
      F32_PO = F_3_2(ETA_PO)
C
      TAU_PO = KQ*F32_PO
C
      TAU_NO = KQ*F32_NO
C
C
C
C
C
C
C                    Calculate internal energy of outside nucleons
c20      BUOUT = (XNUT+XPROT)*( UQ*(TAU_PO+TAU_NO)+
c20     1    EIFLAG*((NOUT**2)*AA+
c20     2   4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM) )/NOUT
      BUOUT = (XNUT+XPROT)*( UQ*(TAU_PO+TAU_NO)+
     1    EIFLAG*PV_E(NPOUT,NNOUT) )/NOUT
C
C
C                                Calc alfa particle internal energy
      BUALFA = 0.25*XALFA*(1.5*T-BALPHA)
C
C                                Set nuclei internal energy to zero
      BUNUC = 0.0
C                                Calculate total baryon internal energy
C                                (per baryon)
      BU = BUOUT+BUALFA+BUNUC
C
C
C                                Calc entropy of outside nucleons
      BSOUT = (XNUT+XPROT)*( (5.0*UQ/(3.0*T))*(TAU_NO+TAU_PO)-
     1   NNOUT*ETA_NO-NPOUT*ETA_PO )/NOUT
C
C                                Calc alpha particle entropy
      BSALFA = 0.25*XALFA*(2.5-MUALFA/T)
C
C                                Set nuclei entropy to zero
      BSNUC = 0.0
C
C                                Calc total baryon entropy (per baryon)
      BS = BSOUT+BSALFA+BSNUC
C
C
C
C                                Calc outside free energy
      BFOUT = BUOUT-T*BSOUT
C                                Calc alpha particle free energy
      BFALFA = BUALFA-T*BSALFA
C                                Set nuclei free energy to zero
      BFNUC = BUNUC-T*BSNUC
C                                Calc total baryon free energy (per nucleon)
      BFTOT = BFOUT+BFALFA+BFNUC
C
C
C
C
C
C                                Calc outside pressure
c20      BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(
c20     1    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)))
      BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
C                                Calc alfa particle pressure
      BPRALF = ALFDNS*T
C
C                                Set nuclei pressure to zero
      BPRNUC = 0.0
C
C                                Calc total baryon pressure
      BPRESS = BPROUT+BPRALF+BPRNUC
C
C
C
C
C
C
C
C
C                           Leptons & Photons
      CALL EL_EOS(T,YE,BRYDNS)
C
C
C
C                           Total pressure and eng/ent per baryon
C
      FBARY = BFTOT+FSUBE
      PBARY = BPRESS+EPRESS
      MUBARY = YE*MUPROT+(1.0-YE)*MUN
      MU_MAT = YE*(MUPROT+MUSUBE)+(1.0-YE)*MUN
C
      FTOT = BFTOT+FSUBE+PF
      UTOT = BU+EU+PU
      STOT = BS+ES+PS
      PTOT = BPRESS+EPRESS+PPRESS
C
C
C
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------
C                Derivatives of thermodynamic variables
C-----------------------------------------------------------------------
C
C
C
cc      GPO = 2.0*FHALFO(ETA_PO)
cc      GNO = 2.0*FHALFO(ETA_NO)
C
C
      GPO = 2.0*FHALF(ETA_PO)
      GNO = 2.0*FHALF(ETA_NO)
C
C
C                 ------------------------------------
C                 !      Derivatives of exterior     !
C                 !      quantities                  !
C                 !      (w.r.t. Temp. and ETA's)    !
C                 !                                  !
C                 ------------------------------------
C
C
C                  Derivatives of exterior potentials
C                  w.r.t. particle densities
c20      DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) )
      DVPODP = EIFLAG*DPVPDP(NPOUT,NNOUT)
c20      DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) )
      DVPODN = EIFLAG*DPVPDN(NPOUT,NNOUT)
c20      DVNODP = DVPODN
      DVNODP = EIFLAG*DPVNDP(NPOUT,NNOUT)
c20      DVNODN = DVPODP
      DVNODN = EIFLAG*DPVNDN(NPOUT,NNOUT)
C
C
C                  Derviatives of exterior chem. pot. w.r.t. ETA's
C                  (at fixed T)
      DMPDEP = T+DVPODP*NPOUT/GPO
      DMPDEN = DVPODN*NNOUT/GNO
      DMNDEP = DVNODP*NPOUT/GPO
      DMNDEN = T+DVNODN*NNOUT/GNO
C
C                  Derivatives of pressure potential w.r.t.
C                  particle densities
c20      DV_DPO = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*DD*(1.0+DD)*(NOUT**DD) )
      DV_DPO = EIFLAG*DPVRDP(NPOUT,NNOUT)
c20      DV_DNO = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NPOUT+CC*DD*(1.0+DD)*(NOUT**DD) )
      DV_DNO = EIFLAG*DPVRDN(NPOUT,NNOUT)
C
C                  Derivatives of pressure potential w.r.t. ETA's
C                  (at fixed T)
      DV_DEP = DV_DPO*NPOUT/GPO
      DV_DEN = DV_DNO*NNOUT/GNO
C
C                  Derivatives of outside pressure w.r.t. ETA's
C                  (at fixed T)
      DPODEP = NPOUT*T+DV_DEP
      DPODEN = NNOUT*T+DV_DEN
C
C                  Derivatives of alpha density w.r.t. ETA's
C                  (at fixed T)
      DNADEP = ALFDNS*(2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP)/T
      DNADEN = ALFDNS*(2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN)/T
C
C
C                  Derivatives of particle densities w.r.t. T
C                  (at fixed ETA's)
      DNPODT = 1.5*NPOUT/T
      DNNODT = 1.5*NNOUT/T
C
C                  Derivatives of exterior chem. pot. w.r.t. T
C                  (at fixed ETA's)
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT
C
C                  Derivative of pressure potential w.r.t. T
C                  (at fixed ETA's)
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT
C
C                  Derivative of outside pressure w.r.t. T
C                  (at fixed ETA's)
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT
C
C                  Derivative of alpha chem. pot. w.r.t. T
C                  (at fixed ETA's)
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT
C
C                  Derivative of alpha particle density w.r.t. T
C                  (at fixed ETA's)
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T
C
C                  Derivative of alpha particle pressure w.r.t. T
C                  (at fixed ETA's)
      DPADT = ALFDNS+T*DNADT
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
C                 ------------------------------------
C                 !      Derivatives of constraint   !
C                 !      and equilibrium equations   !
C                 !      with respect to the five    !
C                 !      compositional variables     !
C                 !      (U,x,n_i,eta_po,eta_no)     !
C                 !      and the three independent   !
C                 !      variables                   !
C                 !      (Baryon density, T, and Ye) !
C                 !                                  !
C                 ------------------------------------
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 1 (Baryon conservation)
C
C
C
      DG1DO1 = -EXALFA*NPOUT/GPO+(V_ALFA*NOUT-4.0)*DNADEP
C
      DG1DO2 = -EXALFA*NNOUT/GNO+(V_ALFA*NOUT-4.0)*DNADEN
C
C
      DG1DL1 = 1.0
C
      DG1DL2 = -EXALFA*(DNNODT+DNPODT)+(V_ALFA*NOUT-4.0)*DNADT
C
      DG1DL3 = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 2 (Charge conservation)
C
C
      DG2DO1 = -EXALFA*NPOUT/GPO+(V_ALFA*NPOUT-2.0)*DNADEP
C
      DG2DO2 = (V_ALFA*NPOUT-2.0)*DNADEN
C
C
      DG2DL1 = YE
C
      DG2DL2 = -EXALFA*DNPODT+(V_ALFA*NPOUT-2.0)*DNADT
C
      DG2DL3 = BRYDNS
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      DET_GT = DG1DO1*DG2DO2-DG1DO2*DG2DO1
C
C
      DEP_DN = (DG1DO2*DG2DL1-DG2DO2*DG1DL1)/DET_GT
      DEN_DN = (DG2DO1*DG1DL1-DG1DO1*DG2DL1)/DET_GT
C
C
      DEP_DT = (DG1DO2*DG2DL2-DG2DO2*DG1DL2)/DET_GT
      DEN_DT = (DG2DO1*DG1DL2-DG1DO1*DG2DL2)/DET_GT
C
C
      DEP_DY = (DG1DO2*DG2DL3-DG2DO2*DG1DL3)/DET_GT
      DEN_DY = (DG2DO1*DG1DL3-DG1DO1*DG2DL3)/DET_GT
C
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of exterior bulk !
C                 !      nucleon internal energy &    !
C                 !      entropy densities and the    !
C                 !      chem. pot.  w.r.t. to eta_p, !
C                 !      ate_n & T. These are used in !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
      S_OUT =(OVR53*UQ/T)*(TAU_NO+TAU_PO)-NNOUT*ETA_NO-NPOUT*ETA_PO
C
c20      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*
c20     1((NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM)
      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*PV_E(NPOUT,NNOUT)
C
C                   Derivative of exterior entropy density w.r.t. T
      DSODT =  OVR53*UQ*(1.5*(TAU_PO+TAU_NO)/(T**2))-
     1     1.5*(NPOUT*ETA_PO+NNOUT*ETA_NO)/T
C
      DEODT = T*DSODT
C
C                    Derivatives of exterior particle densities w.r.t.
C                    Temperature (ETA's fixed)
      DNPODT = 1.5*NPOUT/T
      DNNODT = 1.5*NNOUT/T
C
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT
C
C
      DNPDEP = NPOUT/GPO
      DNNDEN = NNOUT/GNO
C
      DTPDEP = 1.5*T*NPOUT/UQ
      DTNDEN = 1.5*T*NNOUT/UQ
C
      DSODEP = (OVR53*UQ/T)*DTPDEP-NPOUT-ETA_PO*DNPDEP
      DSODEN = (OVR53*UQ/T)*DTNDEN-NNOUT-ETA_NO*DNNDEN
C
C
C                    Exterior particle potentials
c20      VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD )
      VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
c20      VPOUT = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*(1.0+DD)*NOUT**DD+DELTAM)
      VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
      DEODEP = UQ*DTPDEP+VPOUT*DNPDEP
      DEODEN = UQ*DTNDEN+VNOUT*DNNDEN
C
      DMPDEP = T+DVPODP*NPOUT/GPO
      DMPDEN = DVPODN*NNOUT/GNO
      DMNDEP = DVNODP*NPOUT/GPO
      DMNDEN = T+DVNODN*NNOUT/GNO
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of alpha         !
C                 !      particle internal energy &   !
C                 !      entropy densities and the    !
C                 !      chem. pot.  w.r.t. to eta_p, !
C                 !      ate_n & T. These are used in !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
C
C
      S_ALFA = ALFDNS*(2.5-MUALFA/T)
C
C
      E_ALFA = ALFDNS*(1.5*T-BALPHA)
C
C                  Derivative of pressure potential w.r.t. T
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT
C
C                  Derivative of outside pressure w.r.t. T
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT
C
C
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT
C
C                  Derivative of alpha particle density w.r.t. T
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T
C
C
      DSADT = DNADT*(2.5-MUALFA/T)-ALFDNS*DMUADT/T+ALFDNS*MUALFA/T**2
C
      DEADT = DNADT*(1.5*T-BALPHA)+1.5*ALFDNS
C
C
      DV_DEP = DV_DPO*NPOUT/GPO
      DV_DEN = DV_DNO*NNOUT/GNO
C
      DPODEP = OVR23*UQ*DTPDEP+DV_DEP
      DPODEN = OVR23*UQ*DTNDEN+DV_DEN
C
      DMADEP = 2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP
      DMADEN = 2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN
C
      DNADEP = ALFDNS*DMADEP/T
      DNADEN = ALFDNS*DMADEN/T
C
      DSADEP = DNADEP*(2.5-MUALFA/T)-ALFDNS*DMADEP/T
      DSADEN = DNADEN*(2.5-MUALFA/T)-ALFDNS*DMADEN/T
C
      DEADEP = DNADEP*(1.5*T-BALPHA)
      DEADEN = DNADEN*(1.5*T-BALPHA)
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      S_DENS = EXALFA*S_OUT+S_ALFA
C
      E_DENS = EXALFA*E_OUT+E_ALFA
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !      Temperature Derivatives     !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
      DNA_DT = DNADT+DNADEP*DEP_DT+DNADEN*DEN_DT
C
C
      DBSDT = (-V_ALFA*DNA_DT*S_OUT+
     1    EXALFA*(DSODT+DSODEP*DEP_DT+DSODEN*DEN_DT)+
     2    (DSADT+DSADEP*DEP_DT+DSADEN*DEN_DT) )/BRYDNS
C
C~~~~~~~~~~~~~~~~~~
C
      DBUDT = T*DBSDT
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDT = DBUDT-S_DENS/BRYDNS-T*DBSDT
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDT = YE*(DMPODT+DMPDEP*DEP_DT+DMPDEN*DEN_DT)+
     1    (1.0-YE)*(DMNODT+DMNDEP*DEP_DT+DMNDEN*DEN_DT)
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDT = BRYDNS*(DBMUDT-DBFDT)
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !       Density Derivatives        !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
C
      DNA_DN = DNADEP*DEP_DN+DNADEN*DEN_DN
C
C
      DBSDN = (-V_ALFA*DNA_DN*S_OUT+
     1    EXALFA*(DSODEP*DEP_DN+DSODEN*DEN_DN)+
     2   (DSADEP*DEP_DN+DSADEN*DEN_DN) )/BRYDNS-S_DENS/BRYDNS**2
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBUDN = (-V_ALFA*DNA_DN*E_OUT+
     1    EXALFA*(DEODEP*DEP_DN+DEODEN*DEN_DN)+
     2   (DEADEP*DEP_DN+DEADEN*DEN_DN) )/BRYDNS-E_DENS/BRYDNS**2
C
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDN = DBUDN-T*DBSDN
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDN = YE*(DMPDEP*DEP_DN+DMPDEN*DEN_DN)+
     1    (1.0-YE)*(DMNDEP*DEP_DN+DMNDEN*DEN_DN)
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDN = BRYDNS*(DBMUDN-DBFDN)+MUBARY-BFTOT
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !         Ye Derivatives           !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
C
C
C
      DNA_DY = DNADEP*DEP_DY+DNADEN*DEN_DY
C
C
      DBSDY = (-V_ALFA*DNA_DY*S_OUT+
     1    EXALFA*(DSODEP*DEP_DY+DSODEN*DEN_DY)+
     2   (DSADEP*DEP_DY+DSADEN*DEN_DY) )/BRYDNS
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBUDY = (-V_ALFA*DNA_DY*E_OUT+
     1    EXALFA*(DEODEP*DEP_DY+DEODEN*DEN_DY)+
     2   (DEADEP*DEP_DY+DEADEN*DEN_DY) )/BRYDNS
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDY = DBUDY-T*DBSDY
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDY = YE*(DMPDEP*DEP_DY+DMPDEN*DEN_DY)+MUPROT+
     1    (1.0-YE)*(DMNDEP*DEP_DY+DMNDEN*DEN_DY)-MUN
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDY = BRYDNS*(DBMUDY-DBFDY)
C
C
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C-----------------------------------------------------------------------
C                End of derivatives of thermodynamic variables
C-----------------------------------------------------------------------
C
C                  Total derivatives
C                  (Baryons+Electrons+Photons)
C
      DUDT = DBUDT+DEUDT+DPUDT
      DUDN = DBUDN+DEUDN+DPUDN
      DUDY = DBUDY+DEUDY+DPUDY
C
C
      DSDT = DBSDT+DESDT+DPSDT
      DSDN = DBSDN+DESDN+DPSDN
      DSDY = DBSDY+DESDY+DPSDY
C
C
      DPDT = DBPDT+DEPDT+DPPDT
      DPDN = DBPDN+DEPDN+DPPDN
      DPDY = DBPDY+DEPDY+DPPDY
C
C
      DMUDT = DBMUDT+YE*DEMUDT
      DMUDN = DBMUDN+YE*DEMUDN
      DMUDY = DBMUDY+YE*DEMUDY
C
C
C                  Adiabatic index
      GAM_S = BRYDNS*DPDN/PTOT+T*(DPDT**2)/(BRYDNS*PTOT*DUDT)
C
C
      INPVAR(2) = NSUBS
      INPVAR(3) = ETA_PO
      INPVAR(4) = ETA_NO
C
C
C                           Approximate the nuclear density
      NSUBI = NSUBS
C
C                           Use 0.45 as the nuclear proton fraction
      X = 0.45
      A = 4.0D0
C
C                           Save the proton number density for use
C                           as the initial guess on next call
      P_PREV = NPOUT
C
  999 CONTINUE
C
C
      INPVARP(1)     = INPVAR(1)
      INPVARP(2)     = INPVAR(2)
      INPVARP(3)     = INPVAR(3)
      INPVARP(4)     = INPVAR(4)
      YEP            = YE
      BRYDNSP        = BRYDNS
      P_PREVP        = P_PREV
      SSFLAGP        = SSFLAG
C
C
      RETURN
C
      END
Calfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalfaalf

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EL_EOS.FOR
C
C***********************************************************************
C
C    MODULE:       EL_EOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         2/12/91
C
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@SBAST1.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C    PURPOSE:      The elctron and photon equation of state
C
C
C    CALL LINE:    CALL EL_EOS(T,YE,BRYDNS)
C
C    INPUTS:       T = TEMPERATURE
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:      NONE
C
C
C
C    INCLUDE FILES:  EL_EOS.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EL_EOS(T,YE,BRYDNS)
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
C
      IF( ((BRYDNS*YE*1.66D15).GT.5.0D9).OR.(T.GT.5.0D0) ) THEN
C                      Assume that the electrons are relativistic
c        write(*,*) ' rel: '
        CALL EL_REL(T,YE,BRYDNS)
      ELSE
c        write(*,*) ' gl: '
C                      Otherwise call the Gauss-Laguerre version
        CALL EL_GL(T,YE,BRYDNS)
C                      If the stuff is really degenerate call the
C                      relativistic version anyway
        IF((MUSUBE/T).GT.2.0D1) THEN
c        write(*,*) ' rel2: ',musube/t
          CALL EL_REL(T,YE,BRYDNS)
        ENDIF
      ENDIF
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
 999  RETURN
C
C
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EL_GL.FOR
C
C***********************************************************************
C
C    MODULE:       EL_GL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         2/12/91
C                  version 2.0 1/24/93 (includes full fermi integrals)
C
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: swesty@nuclear.physics.sunysb.edu or
C                            fswesty@sbast3.ess.sunysb.edu
C
C    PURPOSE:      The elctron and photon equation of state
C
C
C    CALL LINE:    CALL EL_GL(T,YE,BRYDNS)
C
C    INPUTS:       T = TEMPERATURE
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:      NONE
C
C
C
C    INCLUDE FILES:  EL_EOS.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EL_GL(T,YE,BRYDNS)
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
      DOUBLE PRECISION MASS_E
      PARAMETER(MASS_E=0.51100D0)
C
C                           Pi and 1/3
      DOUBLE PRECISION PI, PI2, OVR3, MOVR3, OVR23
      PARAMETER(PI=3.1415927, PI2=9.8696044)
      PARAMETER(OVR3=0.33333333, MOVR3=-0.33333333, OVR23=0.66666667)
C
C                           2nd Fermi integral
      DOUBLE PRECISION F_2
C                           Positron degeneracy parameter
      DOUBLE PRECISION ELPETA
C                           Multiplicative factors in the G-L
C                           quadrature
      DOUBLE PRECISION TFAC, DTFDT, EFAC, PFAC
C                           Electron eta and beta
      DOUBLE PRECISION ETA_E, BETA, DBTDT, DETDT, DETDNE, D_ETA
C                           Number density and it's derivatives
      DOUBLE PRECISION NE_CHK, NEN, DNENDB, DNENDE
C                           Energy density & derivatives
      DOUBLE PRECISION E_ENG, DEDB, DEDE
C                           Pressure & derivatives
      DOUBLE PRECISION E_PR, DPDB, DPDE
C                           Electron & positron occupation numbers
      DOUBLE PRECISION FE, FP, DFEDB, DFEDE, DFPDB, DFPDE
C                           Integrand factors & derivatives
      DOUBLE PRECISION XN1, DXN1DB, XROOT
      DOUBLE PRECISION XN2, DXN2DB, XN3, DXN3DB
C                           Loop variables
      INTEGER I, J, MAXITR
      PARAMETER(MAXITR=30)
C
C                           Convergence tolerance for N-R iteration
      DOUBLE PRECISION EPSIL
      PARAMETER(EPSIL=1.0D-10)
C
      INTEGER NGAUSS
      PARAMETER(NGAUSS=16)
c      PARAMETER(NGAUSS=36)
      DOUBLE PRECISION WGHT(NGAUSS), X(NGAUSS)
C
      DATA X /.087649410479D00
     *,.46269632892D00,.11410577748D01,.21292836451D01,.34370866339D01
     *,.50780186145D01,.70703385350D01,.94383143364D01,.12214223369D02
     *,.15441527369D02,.19180156857D02,.23515905694D02,.28578729743D02
     *,.34583398702D02,.41940452648D02,.51701160340D02/
      DATA WGHT /.22503631486D00
     *,.52583605276D00,.83196139169D00,.11460992410D01,.14717513170D01
     *,.18131346874D01,.21755175197D01,.25657627502D01,.29932150864D01
     *,.34712344831D01,.40200440864D01,.46725166077D01,.54874206580D01
     *,.65853612333D01,.82763579844D01,.11824277552D02/
C
C23456789012345678901234567890123456789012345678901234567890123456789012
c      DATA X /
c     * 3.961320640860D-02 ,  2.088002856543D-01 ,  5.135107755352D-01 ,
c     * 9.543811530791D-01 ,  1.532249225394D+00 ,  2.248215817854D+00 ,
c     * 3.103661490498D+00 ,  4.100262545699D+00 ,  5.240010108755D+00 ,
c     * 6.525233330468D+00 ,  7.958627508999D+00 ,  9.543288035279D+00 ,
c     * 1.128275128666D+01 ,  1.318104390025D+01 ,  1.524274226656D+01 ,
c     * 1.747304463141D+01 ,  1.987785893308D+01 ,  2.246391051439D+01 ,
c     * 2.523887525730D+01 ,  2.821154567688D+01 ,  3.139204037421D+01 ,
c     * 3.479207144800D+01 ,  3.842529076480D+01 ,  4.230774567208D+01 ,
c     * 4.645849004290D+01 ,  5.090042150332D+01 ,  5.566145791468D+01 ,
c     * 6.077624068453D+01 ,  6.628869068025D+01 ,  7.225601475446D+01 ,
c     * 7.875533816278D+01 ,  8.589548142950D+01 ,  9.383992978484D+01 ,
c     * 1.028580101471D+02 ,  1.134651354897D+02 ,  1.269880151966D+02 /
C
c      DATA WGHT /
c     * 1.016694751217D-01 ,  2.368502055548D-01 ,  3.726715786250D-01 ,
c     * 5.092083225550D-01 ,  6.467110475491D-01 ,  7.854519690232D-01 ,
c     * 9.257184440339D-01 ,  1.067815135684D+00 ,  1.212067715646D+00 ,
c     * 1.358827381029D+00 ,  1.508476196330D+00 ,  1.661433435563D+00 ,
c     * 1.818163187438D+00 ,  1.979183577131D+00 ,  2.145078075547D+00 ,
c     * 2.316509528323D+00 ,  2.494237764181D+00 ,  2.679141968116D+00 ,
c     * 2.872249479966D+00 ,  3.074773384330D+00 ,  3.288162326323D+00 ,
c     * 3.514167642788D+00 ,  3.754935526960D+00 ,  4.013136236755D+00 ,
c     * 4.292149591579D+00 ,  4.596338646060D+00 ,  4.931466465190D+00 ,
c     * 5.305354988656D+00 ,  5.728974249695D+00 ,  6.218344086884D+00 ,
c     * 6.798089006814D+00 ,  7.508695084046D+00 ,  8.423188762282D+00 ,
c     * 9.692562072805D+00 ,  1.170765394670D+01 ,  1.607328906573D+01 /
C
C
C-----------------------------------------------------------------------
C
C
C
C
C                    Electron number density
      NSUBE = BRYDNS*YE
C
C                    Coefficants for chemical potential
C                    and thermodynamics quantities
      QSUBE = 1.0D0/( 3.0D0*(PI**2)*((HBAR*C)**3) )
C
      ACOEF = 0.5D0*NSUBE/QSUBE
C
      BCOEF = (ACOEF**2+((PI**6)*T**6)/27.0D0)**0.5D0
C
      DBDT = (PI**6)*(T**5)/(9.0D0*BCOEF)
C
      CCOEF = (ACOEF+BCOEF)**OVR3
C
C
C                    Electron chemical potential
      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF
C
C
C
C                            Initial guess at electron degeneracy
C                            parameter
      ETA_E = MUSUBE/T
C
      BETA = MASS_E/T
      DBTDT = -BETA/T
C                            Multiplicative factor in Fermi integrals
      TFAC = (T**3)/((PI**2)*(HBAR*C)**3)
      DTFDT = 3.0D0*TFAC/T
      PFAC = TFAC*T*OVR3
      EFAC = TFAC*T/BRYDNS
C
C                            Loop until the N-R iteration for the
C                            chemical potential converges
      DO 20 J=1,MAXITR,1
C
C                            Zero the electron number density integral
        NEN = 0.0D0
        DNENDE = 0.0D0
        DNENDB = 0.0D0
C
C                            Zero the electron energy density integral
        E_ENG = 0.0D0
        DEDE = 0.0D0
        DEDB = 0.0D0
C
C                            Zero the electron pressure integral
        E_PR = 0.0D0
        DPDE = 0.0D0
        DPDB = 0.0D0
C
C                            Zero the positron number density
        NEPLUS = 0.0D0
C
C       ----------------------------------------------------------------
C       |                    Do the Gaussian Quadrature                |
        DO 10 I=1,NGAUSS,1
C       |                    Electron occupation numbers               |
          FE = 1.0D0/(DEXP(X(I)+BETA-ETA_E)+1.0D0)
C       |                    Positron occupation numbers               |
          FP = 1.0D0/(DEXP(X(I)+BETA+ETA_E)+1.0D0)
C       |                    Derivatives w.r.t. beta                   |
          DFEDB = -FE*(1.0D0-FE)
          DFPDB = -FP*(1.0D0-FP)
C       |                    Derivatives w.r.t. Eta                    |
          DFEDE = FE*(1.0D0-FE)
          DFPDE = -FP*(1.0D0-FP)
C       |                                                              |
C       |                    Number density integral & derivatives     |
          XROOT = DSQRT(X(I)*(X(I)+2.0D0*BETA))
          XN1 = XROOT*(X(I)+BETA)
          DXN1DB = XROOT+X(I)*(X(I)+BETA)/XROOT
          NEN = NEN+WGHT(I)*XN1*(FE-FP)
          DNENDB = DNENDB+WGHT(I)*(DXN1DB*(FE-FP)+XN1*(DFEDB-DFPDB))
          DNENDE = DNENDE+WGHT(I)*XN1*(DFEDE-DFPDE)
C       |                    Sum for positron # integral               |
          NEPLUS = NEPLUS+WGHT(I)*XN1*FP
C       |                                                              |
C       |                                                              |
C       |                    Energy integral & derivatives             |
          XN2 = XROOT*(X(I)+BETA)**2
          DXN2DB = 2.0D0*XROOT*(X(I)+BETA)+
     1          X(I)*((X(I)+BETA)**2)/XROOT
          E_ENG = E_ENG+WGHT(I)*XN2*(FE+FP)
          DEDB = DEDB+WGHT(I)*(DXN2DB*(FP+FE)+XN2*(DFEDB+DFPDB))
          DEDE = DEDE+WGHT(I)*XN2*(DFEDE+DFPDE)
C       |                                                              |
C       |                                                              |
C       |                    Pressure integral & derivatives           |
          XN3 = XROOT**3
          DXN3DB = 3.0D0*XROOT*X(I)
          E_PR = E_PR+WGHT(I)*XN3*(FE+FP)
          DPDB = DPDB+WGHT(I)*(DXN3DB*(FE+FP)+XN3*(DFEDB+DFPDB))
          DPDE = DPDE+WGHT(I)*XN3*(DFEDE+DFPDE)
C       |                                                              |
C       |                                                              |
 10     CONTINUE
C       |                                                              |
C       ----------------------------------------------------------------
C                            Multiply by the temperature factor to get
C                            number density
        NE_CHK = TFAC*NEN
C                            Multiply by the temperature factor to get
C                            number density
        NEPLUS = TFAC*NEPLUS
C
C                            Calculate the new change in ETA
        D_ETA = -(NE_CHK-NSUBE)/(DNENDE*TFAC)
C
C                   If we've met the convergence criterion...
        IF(DABS(D_ETA).LT.(EPSIL*ETA_E)) THEN
C                          Then break out the chemical potential
C                          loop
          GOTO 30
        ELSE
C                          Otherwise update the chemical potential
          ETA_E = ETA_E+D_ETA
        ENDIF
C
 20   CONTINUE
C
C                                      If we reached this point the
C                                      N-R iteration didn't converge
      WRITE(*,*) ' EL_EOS: N-R iteration didnt converge!, 
     * T,YE,BRYDNS=',T,YE,BRYDNS
      STOP
C
C               If we reached this point then the N-R iteration has
C               converged.
 30   CONTINUE
C
C                     Is the result consistent with electron number
C                     density?
      IF(DABS((NE_CHK-NSUBE)/NSUBE).GT.1.0D-5) THEN
        WRITE(*,*) ' EL_EOS: Gaussian quadrature converged badly!'
        STOP
      ENDIF
C
C             Calculate thermodynamic quantities...
C
C                          Electron chemical potential
      MUSUBE = T*ETA_E
C                          Electron pressure
      EPRESS = PFAC*E_PR
C                          Electron internal energy per baryon
      EU = EFAC*E_ENG
C                          Electron free energy per baryon
      FSUBE = YE*MUSUBE-EPRESS/BRYDNS
C                          Electron entropy per baryon
      ES = (EU-FSUBE)/T
C
C                          Derivative of the electron eta w.r.t. T
      DETDT = -(DTFDT*NEN+TFAC*DNENDB*DBTDT)/(TFAC*DNENDE)
C                          Derivative of the electron eta w.r.t. NSUBE
      DETDNE = 1.0D0/(TFAC*DNENDE)
C
C                    Derivatives of chem. potential w.r.t. T,
C                    BRYDNS, YE
      DEMUDT = T*DETDT+MUSUBE/T
      DEMUDN = YE*T*DETDNE
      DEMUDY = BRYDNS*T*DETDNE
C
C
C                    Derivatives of pressure w.r.t. BRYDNS,YE,T
      DEPDN = YE*PFAC*DPDE*DETDNE
      DEPDY = BRYDNS*PFAC*DPDE*DETDNE
      DEPDT = (4.0D0*PFAC/T)*E_PR+PFAC*(DPDE*DETDT+DPDB*DBTDT)
C
C                    Derivatives of internal energy w.r.t.
C                    T,BRYDNS,YE
      DEUDT = (4.0D0*EFAC/T)*E_ENG+EFAC*(DEDB*DBTDT+DEDE*DETDT)
      DEUDN = (-1.0D0*EFAC/BRYDNS)*E_ENG+YE*EFAC*DEDE*DETDNE
      DEUDY = BRYDNS*EFAC*DEDE*DETDNE
C
C                    Derivatives of entropy w.r.t. T,BRYDNS,YE
      DESDT = -ES/T+(DEUDT-YE*DEMUDT+DEPDT/BRYDNS)/T
      DESDN = (DEUDN-YE*DEMUDN+DEPDN/BRYDNS-EPRESS/(BRYDNS**2))/T
      DESDY = (DEUDY-MUSUBE-YE*DEMUDY+DEPDY/BRYDNS)/T
C
C
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
C23456789012345678901234567890123456789012345678901234567890123456789012
C
 999  RETURN
C
C
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       F_2
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         12/16/91
C                  Version 2: 1/24/93
C
C    CALL LINE:    F_2(Y)      (2nd Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       2nd Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION F_2(Y)
      IMPLICIT NONE
      DOUBLE PRECISION Y, YEXP
      IF(Y.GT.3.0D0) THEN
        WRITE(*,*) ' F_2(Y) FAILS FOR Y .GT. 3; Y =',Y
        STOP
      ENDIF
C
C                       Note: This approximation is based on the
C                       Bludman & Van Riper approximation (see
C                       Ap. J. Vol. 212 page 866-867 (1977))
C                       equation (3.6)
C
      IF(Y.LT.-1.0D0) THEN
        YEXP = EXP(Y)
        F_2 = 2.0*YEXP*(1.0-0.125*YEXP+0.037037*(YEXP**2))
      ELSE
        F_2 = 1.803D0+1.645D0*Y+0.6931D0*(Y**2)+0.1666667D0*(Y**3)+
     1        2.0833333D-2*(Y**4)-3.4722D-4*(Y**6)
      ENDIF
 999  RETURN
      END

!///////////////////////////////////////////////////////////////////////

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

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       EOS_M4C
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         3/3/92  Model 4C modifications completed
C                  12/15/90 Modified from model 4A to include the
C                  phase boundary cutoffs and Maxwell construction
C                  boundaries.
C                  7/13/90 Modified from model 1-d to include Maxwell
C                  construction
C                  5/25/90  MODEL 1D
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C
C    CALL LINE:    CALL EOS_M4C(INPVAR,YE,BRYDNS,IFLAG,EOSFLG,FFLAG,
C                  XPREV,P_PREV)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C                  IFLAG = 1 --> INPVAR IS TEMPERATURE
C                          2 --> INPVAR IS INTERNAL ENERGY (NOT IMPLEM)
C                          3 --> INPVAR IS ENTROPY (NOT IMPLEMENTED)
C                          (IFLAG=1 is now assumed at this level)
C                  FFLAG = "FORCING FLAG"  0 --> NO FORCING
C                                          1 --> FORCE A PARTICULAR
C                                                SCHEME TO BE USED
C
C
C    OUTPUTS:      EOSFLG = 1 --> Not implemented in model 4B
C                           2 --> GENERAL EOS
C                           3 --> BULK EOS (includes alpha's)
C                  XPREV = PREVIOUS VALUE OF X (MUST BE SUPPLIED ON
C                          FIRST CALL)
C                  P_PREV = PREVIOUS VALUE OF PROTON DENSITY (MUST BE
C                          SUPPLIED ON FIRST CALL)
C
C
C
C
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EOS_M4C(INPVARP,YEP,BRYDNSP,IFLAGP,EOSFLGP,FFLAGP,
     1                   SSFLAGP,XPREVP,P_PREVP)
C
      USE eos_m4c_module
C
      IMPLICIT NONE
C
      DOUBLE PRECISION OUTVAR(4)
C
      INTEGER IFLAGP,EOSFLGP,FFLAGP,SSFLAGP
      DOUBLE PRECISION INPVARP(4),YEP,BRYDNSP,XPREVP,P_PREVP
C
C
      INPVAR(1)      = INPVARP(1)
      INPVAR(2)      = INPVARP(2)
      INPVAR(3)      = INPVARP(3)
      INPVAR(4)      = INPVARP(4)
      YE             = YEP
      BRYDNS         = BRYDNSP
      IFLAG          = IFLAGP
      EOSFLG         = EOSFLGP
      FFLAG          = FFLAGP
      SSFLAG         = SSFLAGP
      XPREV          = XPREVP
      P_PREV         = P_PREVP
C
C
C                       This include file contains all variable
C                       declarations.  NOTE:: no implicit typing
C                       scheme is followed in this code; if you
C                       have any doubt as to a variables type CHECK
C                       IT!!!!.  Also note that ALL variables are
C                       declared explicitly.
C
C
C                         Set the "switch" flag to zero
      SWTFLG = 0
C
C                         Set T equal to the input variable (the entropy
C                         and internal energy options should go through
C                         INVEOS untill further notice)
      T = INPVAR(1)
C
C
C                         If the "forcing" flag is set then skip
C                         the EOS determination logic and go straight
C                         to the EOS determined by EOSFLG
      IF(FFLAG.EQ.1) THEN
        GOTO 10
      ELSE
C                         Otherwise let the EOS logic module determine
C                         the correct EOS to use
        CALL EOSLOG(INPVAR,YE,BRYDNS,EOSFLG)
      ENDIF
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                        Try NUCEOS first and if not successfull
C                        then try bulk EOS
 10   CONTINUE
      IF(EOSFLG.EQ.1) THEN
C
        CALL NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C                    If the nuclear EOS failed and the reset flag is set
C                    then reset the initial guesses and try again
        IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN
          CALL RESET(INPVAR,YE,BRYDNS,OUTVAR)
          OUTVAR(1) = INPVAR(1)
          CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C
C                    Make a last ditch effort at convergence
          IF(SSFLAG.NE.1) THEN
            OUTVAR(2) = 0.155
            OUTVAR(3) = -15.0
            OUTVAR(4) = -20.0
            CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
          ENDIF
C
        ENDIF
C
C
C
        IF((XH.GT.HEAVCT).AND.(SSFLAG.EQ.1)) THEN
C                    Set EOS flag to full scheme
          EOSFLG = 2
C
C                    Else if fraction of nuclei is less than the minimum
C                    or if NUCEOS was unsuccessful use the no nuclei EOS
        ELSE
          IF(FFLAG.NE.1) THEN
C
            CALL ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG)
C
            IF((SSFLAG.NE.1).AND.(FFLAG.EQ.1)) THEN
              EOSFLG = 1
              WRITE(*,*) 'A2 failed at try = ',T,BRYDNS,YE
              GOTO 999
            ENDIF
C
C                    Set nuclei to bulk EOS
            EOSFLG = 3
C                    Save value of proton fraction
            P_PREV = YE*BRYDNS
C
            GOTO 999
C
          ELSE
            IF(NF_FLG.EQ.1)
     1          WRITE(*,*) 'NUC failed at t,rho = ',t,brydns
            GOTO 999
          ENDIF
        ENDIF
C
      ENDIF
C
C
C          End of NUCEOS--BULK EOS calculations
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C
C
C
C
C
C
C
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                            CALCULATE FULL EOS (INCLUDING NUCLEI)
      IF(EOSFLG.EQ.2) THEN
C
C                    Call the nuclear EOS
        CALL NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C
C                    If the nuclear EOS failed and the reset flag is set
C                    then reset the initial guesses and try again
        IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN
cccc          WRITE(*,*) ' EOS_M4C:: r.i.gs.'
          CALL RESET(INPVAR,YE,BRYDNS,OUTVAR)
          OUTVAR(1) = INPVAR(1)
          CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C
C                    Make a last ditch effort at convergence
          IF(SSFLAG.NE.1) THEN
            OUTVAR(2) = 0.155
            OUTVAR(3) = -15.0
            OUTVAR(4) = -20.0
            CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
          ENDIF
C
C
C
          IF(SSFLAG.NE.1) THEN
cccc            WRITE(*,*) '     r.i.gs. failure @ try: ',inpvar
            GOTO 999
          ELSE
            INPVAR(2) = OUTVAR(2)
            INPVAR(3) = OUTVAR(3)
            INPVAR(4) = OUTVAR(4)
          ENDIF
C                    Otherwise quit and return
        ELSEIF((SSFLAG.NE.1).AND.(FFLAG.EQ.1)) THEN
          GOTO 999
        ENDIF
C
C
C
C
C                    If fraction of heavies is greater than the minimum
C                    parameter, then this EOS is OK
        IF((XH.GT.HEAVCT).AND.(SSFLAG.EQ.1)) THEN
C                    Set EOS flag to full scheme
          EOSFLG = 2
C
C                    Else if fraction of nuclei is less than the minimum
C                    or if NUCEOS was unsuccessful use the no nuclei EOS
        ELSE
C                    If the forcing flag is not set
          IF(FFLAG.NE.1) THEN
C                    Set nuclei to no nuclei EOS
            EOSFLG = 3
C                    Set flag to indicate switch is being made
            SWTFLG = 1
C
            WRITE(*,*) ' NUCEOS failed at try =',t,brydns,ye
            WRITE(*,*) ' where it shouldnt have; Bulk EOS was used'
            WRITE(*,*) ' IV = ',INPVAR
            WRITE(*,*) ' '
C
C                    Branch to bulk EOS
            GOTO 50
C
C                    Otherwise since forcing flag is set then declare
C                    a failure and return
          ELSE
C                      If the failure message flag is set then announce
C                      the failure
            IF(NF_FLG.EQ.1)
     1          WRITE(*,*) 'NUC failed at t,r = ',t,brydns
            GOTO 999
          ENDIF
        ENDIF
C
      ENDIF
C                              END OF FULL EOS CALULATIONS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C
C
C
C
C
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                              CALCULATE BULK EOS
 50   CONTINUE
      IF(EOSFLG.EQ.3) THEN
C
        CALL ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG)
C
        IF((SSFLAG.EQ.0).AND.(FFLAG.EQ.1).AND.(NF_FLG.EQ.1)) THEN
          WRITE(*,*) 'A1 failed at t,rho = ',t,brydns
          GOTO 999
        ENDIF
C                           If this EOS was used as a result of the
C                           nuclear EOS failing then set the
C                           success flag to indicate a warning
        IF(SWTFLG.EQ.1) THEN
          SSFLAG = 2
        ENDIF
C
C                           Save the value of the proton fraction
        P_PREV = YE*BRYDNS
C
        GOTO 999
C
      ENDIF
C                END OF BULK EOS CALCULATIONS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C
C
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                              CALCULATE VIA MAXWELL CONSTRUCTION
      IF(EOSFLG.EQ.4) THEN
C
        CALL MAXWEL(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG)
C
C                 Save the value of the proton fraction
        P_PREV = YE*BRYDNS
C
C                 If Maxwell EOS failed then announce the failure
        IF(SSFLAG.NE.1) THEN
          WRITE(*,*) ' MAXWEL failed at try = '
          WRITE(*,*) T,BRYDNS,YE
        ENDIF
C
          GOTO 999
C
      ENDIF
C                END OF MAXWELL CONSTRUCTION CALCULATIONS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
  999 CONTINUE
C
C
      INPVARP(1)     = INPVAR(1)
      INPVARP(2)     = INPVAR(2)
      INPVARP(3)     = INPVAR(3)
      INPVARP(4)     = INPVAR(4)
      YEP            = YE
      BRYDNSP        = BRYDNS
      IFLAGP         = IFLAG
      EOSFLGP        = EOSFLG
      FFLAGP         = FFLAG
      SSFLAGP        = SSFLAG
      XPREVP         = XPREV
      P_PREVP        = P_PREV
C
C
      RETURN
C
      END
C

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EOSLOG.FOR
C
C***********************************************************************
C
C    MODULE:       EOSLOG
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         12/15/90
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C
C    CALL LINE:    CALL EOSLOG(INPVAR,YE,BRYDNS,EOSFLG)
C
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C
C
C    OUTPUTS:      EOSFLG = 1 --> Not implemented in model 4B
C                           2 --> GENERAL EOS
C                           3 --> BULK EOS (includes alpha's)
C
C
C
C
C    INCLUDE FILES:  EOS_M4C.INC, MAXWEL.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE EOSLOG(INPVARP,YEP,BRYDNSP,EOSFLGP)
C
C
      USE eos_m4c_module
      USE maxwel_module
C
C
      IMPLICIT NONE
C
C                       This include file contains all variable
C                       declarations.  NOTE:: no implicit typing
C                       scheme is followed in this code; if you
C                       have any doubt as to a variables type CHECK
C                       IT!!!!.  Also note that ALL variables are
C                       declared explicitly.
C
C
      DOUBLE PRECISION NLOW, NHI, N_CUT, TEMP_1, TEMP_2, T_BNDY
C
      DOUBLE PRECISION LMM, LMP, LPM, LPP
      DOUBLE PRECISION DNDY1, DNDY2
C
      INTEGER EOSFLGP
      DOUBLE PRECISION INPVARP(4),YEP,BRYDNSP
C
C
C
      INPVAR(1)      = INPVARP(1)
      INPVAR(2)      = INPVARP(2)
      INPVAR(3)      = INPVARP(3)
      INPVAR(4)      = INPVARP(4)
      YE             = YEP
      BRYDNS         = BRYDNSP
      EOSFLG         = EOSFLGP
C
C
 10   CONTINUE
C
C
C                         Set T equal to the input variable (any calls
C                         with entropy or internal energy should go
C                         the the EOS_M4B subroutine)
C
      T = INPVAR(1)
C
C-----------------------------------------------------------------------
C         code to figure out the boundaries from the tables
C-----------------------------------------------------------------------
C
C
C
C
      IF(YE.GT.Y_HI) THEN
C                         Ye is too large for EOS
C
        WRITE(*,*) ' EOSLOG:: Cant do Ye = ',YE, 'at this time'
        WRITE(*,*) ' EOSLOG:: assuming YE =',Y_HI,' instead'
        YE = Y_HI-1.0D-6
        GOTO 10
C
      ELSEIF(YE.GE.Y_LOW) THEN
C                         Calculate high and low boundary densities
C                         for the Maxwell construction
C
C----------------------------------------------------------
C           Calc Ye index
C----------------------------------------------------------
C
        YFRAC = (YE-Y_LOW)/(Y_HI-Y_LOW)
        J_MXWL = INT(YFRAC*(NUMYE-1))+1
        DELT_Y = (Y_HI-Y_LOW)/DBLE(NUMYE-1)
C
        YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y
        YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y
C
C
        IF((YE.GE.YMINUS).AND.(YE.LE.YPLUS)) THEN
          J_BD = J_MXWL
          J_BNDY = J_MXWL
        ELSEIF(YE.GT.YPLUS) THEN
          J_MXWL = J_MXWL+1
          J_BD = J_MXWL
          J_BNDY = J_MXWL
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y
        ELSE
          J_MXWL = J_MXWL-1
          J_BD = J_MXWL
          J_BNDY = J_MXWL
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y
        ENDIF
C
C
        IF(J_MXWL.GT.(NUMYE-1)) THEN
          J_MXWL = NUMYE-1
          J_BD = J_MXWL
          J_BNDY = J_MXWL
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y
        ENDIF
C
C
        YINTRP = (YE-YMINUS)/(YPLUS-YMINUS)
C
C----------------------------------------------------------
C           Calc T index
C----------------------------------------------------------
C
C
        TFRAC = (T-T_LOW)/(T_HI-T_LOW)
        I_MXWL = INT(TFRAC*(NUMTMP-1))+1
        DELT_T = (T_HI-T_LOW)/DBLE(NUMTMP-1)
C
        TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
        TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
C
C
        IF((T.GT.TMINUS).AND.(T.LE.TPLUS)) THEN
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
        ELSEIF(T.GT.TPLUS) THEN
          I_MXWL = I_MXWL+1
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
        ELSE
          I_MXWL = I_MXWL-1
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
        ENDIF
C
C
        IF(I_MXWL.GT.(NUMTMP-1)) THEN
          I_MXWL = NUMTMP-1
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T
        ENDIF
C
C
        TINTRP = (T-TMINUS)/(TPLUS-TMINUS)
C
C
C
C
C                Find the temperature and density at the top of the
C                Maxwel construction
C
CC      T_MXWL = YINTRP*(T_H(J_MXWL+1)-T_H(J_MXWL))+T_H(J_MXWL)
CC      D_MXWL = YINTRP*(D_H(J_MXWL+1)-D_H(J_MXWL))+D_H(J_MXWL)
        T_MXWL = DMIN1(T_H(J_MXWL+1),T_H(J_MXWL))
        IF(T_H(J_MXWL+1).GT.T_H(J_MXWL)) THEN
          D_MXWL = D_H(J_MXWL)
        ELSE
          D_MXWL = D_H(J_MXWL+1)
        ENDIF
C
C
C
C--------------------------------------------------------------------
C            Interpolate to get Maxwell construction densities
C--------------------------------------------------------------------
C
C
C
        DNS_1 = YINTRP*(BRYLOW(I_MXWL,J_MXWL+1)-BRYLOW(I_MXWL,J_MXWL))+
     1               BRYLOW(I_MXWL,J_MXWL)
        DNS_2 = YINTRP*
     1        (BRYLOW(I_MXWL+1,J_MXWL+1)-BRYLOW(I_MXWL+1,J_MXWL))+
     2               BRYLOW(I_MXWL+1,J_MXWL)
C
        LOWDNS = TINTRP*(DNS_2-DNS_1)+DNS_1
C
C                Derivative of lower density w.r.t. T
        DNL_DT = (DNS_2-DNS_1)/DELT_T
C
        DNDY1 = (BRYLOW(I_MXWL,J_MXWL+1)-BRYLOW(I_MXWL,J_MXWL))/DELT_Y
        DNDY2 = (BRYLOW(I_MXWL+1,J_MXWL+1)-
     1      BRYLOW(I_MXWL+1,J_MXWL))/DELT_Y
        DNL_DY = TINTRP*(DNDY2-DNDY1)+DNDY1
C
C
C
C
        IF(YE.GT.Y_CUT) THEN
C
          DNS_1 = YINTRP*
     1        (BRYHI(I_MXWL,J_MXWL+1)-BRYHI(I_MXWL,J_MXWL))+
     2        BRYHI(I_MXWL,J_MXWL)
          DNS_2 = YINTRP*
     1        (BRYHI(I_MXWL+1,J_MXWL+1)-BRYHI(I_MXWL+1,J_MXWL))+
     2               BRYHI(I_MXWL+1,J_MXWL)
C
          HIDNS = TINTRP*(DNS_2-DNS_1)+DNS_1
C
C                Derivative of higher density w.r.t. T
          DNH_DT = (DNS_2-DNS_1)/DELT_T
C
C
        DNDY1 = (BRYHI(I_MXWL,J_MXWL+1)-
     1      BRYHI(I_MXWL,J_MXWL))/DELT_Y
        DNDY2 = (BRYHI(I_MXWL+1,J_MXWL+1)-
     1      BRYHI(I_MXWL+1,J_MXWL))/DELT_Y
        DNH_DY = TINTRP*(DNDY2-DNDY1)+DNDY1
C
C
        ELSE
          HIDNS = LOWDNS
        ENDIF
C
C
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C
C                       Ye is too low
      ELSE
        WRITE(*,*) ' EOSLOG:: Cant do Ye = ',YE, 'at this time'
        WRITE(*,*) ' EOSLOG:: assuming YE =',Y_LOW,' instead'
        YE = Y_LOW+1.0D-6
        GOTO 10
      ENDIF
C
C
C
C
C
      DLTLN1 = (LNCUT-LNLOW)/DBLE(NUMLOW-1)
      DLTLN2 = (LNHI-LNCUT)/DBLE(NUMHI-1)
C
C
      NLOW = 10.0**LNLOW
      NHI = 10.0**LNHI
      N_CUT = 10.0**LNCUT
      LOGBRY = DLOG10(BRYDNS)
      LOGBCH = LOGBRY
C
C
C
C
C----------------------------------------------------------
C           Calc T index
C----------------------------------------------------------
C
C
      IF(LOGBRY.GE.LNHI) THEN
        I_BD = NBPNTS
        I_BNDY = NBPNTS
        T_BNDY = YINTRP*
     1           (LBOUND(I_BNDY,J_BNDY+1)-LBOUND(I_BNDY,J_BNDY))+
     2            LBOUND(I_BNDY,J_BNDY)
        GOTO 70
      ELSEIF((LOGBRY.LT.LNHI).AND.(LOGBRY.GT.LNCUT)) THEN
C
        I_BD = INT((LOGBRY-LNCUT)/DLTLN2)+NUMLOW
        LNMINS = LNCUT+DBLE(I_BD-NUMLOW)*DLTLN2
        LNPLUS = LNCUT+DBLE(I_BD-NUMLOW+1)*DLTLN2
        IF((LOGBCH.LE.LNPLUS).AND.(LOGBCH.GE.LNMINS)) THEN
          I_BNDY = I_BD
        ELSEIF(LOGBCH.GT.LNPLUS) THEN
          I_BD = I_BD+1
          I_BNDY = I_BD
          LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2
          LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2
        ELSE
          I_BD = I_BD-1
          I_BNDY = I_BD
          LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2
          LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2
        ENDIF
C
      ELSEIF((LOGBRY.LE.LNCUT).AND.(LOGBRY.GT.LNLOW)) THEN
C
        I_BD = INT((LOGBRY-LNLOW)/DLTLN1)+1
        LNMINS = LNLOW+DBLE(I_BD-1)*DLTLN1
        LNPLUS = LNLOW+DBLE(I_BD)*DLTLN1
        IF((LOGBCH.LE.LNPLUS).AND.(LOGBCH.GE.LNMINS)) THEN
          I_BNDY = I_BD
        ELSEIF(LOGBCH.GT.LNPLUS) THEN
          I_BD = I_BD+1
          I_BNDY = I_BD
          LNMINS = LNLOW+DBLE(I_BNDY-1)*DLTLN1
          LNPLUS = LNLOW+DBLE(I_BNDY)*DLTLN1
        ELSE
          I_BD = I_BD-1
          I_BNDY = I_BD
          LNMINS = LNLOW+DBLE(I_BNDY-1)*DLTLN1
          LNPLUS = LNLOW+DBLE(I_BNDY)*DLTLN1
        ENDIF
C
      ENDIF
C
      IF(I_BNDY.GT.(NBPNTS-1)) THEN
        I_BD = NBPNTS-1
        I_BNDY = I_BD
        LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2
        LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2
      ENDIF
C
C
C
      LMM = LBOUND(I_BNDY,J_BNDY)
      LPM = LBOUND(I_BNDY+1,J_BNDY)
      LMP = LBOUND(I_BNDY,J_BNDY+1)
      LPP = LBOUND(I_BNDY+1,J_BNDY+1)
C
      LNFRAC = (LOGBCH-LNMINS)/(LNPLUS-LNMINS)
C
C                Interpolate in Ye first
C
      TEMP_1 = YINTRP*
     1           (LBOUND(I_BNDY,J_BNDY+1)-LBOUND(I_BNDY,J_BNDY))+
     2            LBOUND(I_BNDY,J_BNDY)
      TEMP_2 = YINTRP*
     1        (LBOUND(I_BNDY+1,J_BNDY+1)-LBOUND(I_BNDY+1,J_BNDY))+
     2               LBOUND(I_BNDY+1,J_BNDY)
C
C                Interpolate in density between the two Ye
C                interpolated values
C
      T_BNDY = LNFRAC*(TEMP_2-TEMP_1)+TEMP_1
C
C
C----------------------------------------------------------
C----------------------------------------------------------
C
 70   CONTINUE
C
      TCHK_B = 1.01*T_BNDY
      TCHK_N = 0.95*T_BNDY
C
      IF((LMM.GE.LPM).OR.(LMP.GT.LPP)) THEN
        TCHK_N = DMAX1(0.0D0,DMIN1(0.95*TCHK_N,T_BNDY-3.0))
      ENDIF
C
C-----------------------------------------------------------------------
C               EOS Logic
C-----------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C                     If T is below the maximum for maxwel construction
      IF(T.LT.T_MXWL) THEN
C                       If rho is greater than the upper max. con.
C                       density the use the bulk EOS
        IF(BRYDNS.GT.HIDNS) THEN
          EOSFLG = 3
C                       Else if rho is greater than the lower max. con.
C                       density then
        ELSEIF(BRYDNS.GT.LOWDNS) THEN
C                         If Ye is large enough to have a signifigant
C                         max con then use the maxwell con. EOS
          IF(YE.GT.Y_CUT) THEN
            EOSFLG = 4
C                         Otherwise use the bulk EOS
          ELSE
            EOSFLG = 3
          ENDIF
C
C                       If density is greater than the minimum
C                       Maxwell con. density, then we know that we are
C                       in the Nuclear EOS density
        ELSEIF(BRYDNS.GT.D_MXWL) THEN
          EOSFLG = 2
C
C
C                       Otherwise check the Boundary table
        ELSE
C
C                         If T is well below the phase boundary curve
C                         then use the nuclear EOS
          IF(T.LT.TCHK_N) THEN
            EOSFLG = 2
C                         Otherwise if T is near the boundary, first
C                         try the nuclear EOS and if not successfull
C                         then use the bulk EOS
          ELSEIF(T.LT.TCHK_B) THEN
            EOSFLG = 1
          ELSE
C                         Otherwise T is well above the boundary so
C                         use the bulk EOS
            EOSFLG = 3
          ENDIF
        ENDIF
C
C                     Otherwise T is above the maximum for a maxwell
C                     construction
      ELSE
C                       If density is greater than that at the top of
C                       the maxwell construction then use the bulk EOS
        IF(BRYDNS.GT.D_MXWL) THEN
          EOSFLG = 3
C
C                       Otherwise density is below the maxwell con.
        ELSE
C
C                         If T is well below the phase boundary curve
C                         then use the nuclear EOS
          IF(T.LT.TCHK_N) THEN
            EOSFLG = 2
C
C                         Otherwise if T is near the phase boundary
C                         curve then try the nuclear EOS and if not
C                         successfull then use the bulk EOS
          ELSEIF(T.LT.TCHK_B) THEN
            EOSFLG = 1
C
C                         Otherwise T is well above the phase boundary
C                         curve so use the bulk EOS
          ELSE
            EOSFLG = 3
          ENDIF
        ENDIF
      ENDIF
C
C
C-----------------------------------------------------------------------
C                         Done with EOS logic so return EOSFLG
C-----------------------------------------------------------------------
C
C
  999 CONTINUE
C
      INPVARP(1)     = INPVAR(1)
      INPVARP(2)     = INPVAR(2)
      INPVARP(3)     = INPVAR(3)
      INPVARP(4)     = INPVAR(4)
      YEP            = YE
      BRYDNSP        = BRYDNS
      EOSFLGP        = EOSFLG
C
      RETURN
C
C
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       F_1_2
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         11/29/89
C
C    CALL LINE:    F_1_2(Y)      (1/2th Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       1/2th Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION F_1_2(Y)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(n=201)
      DIMENSION A(7),eta(n),f32(n),f32a(n),f12(n),fr(n),f12a(n),fra(n),
     > fia(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
      DATA A,th,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0
     1,6.41500299D-02,.4D0,1.32934039D0,.33333333333d0,1,n/
      IF(y .gt. 30.) goto 10
      if(y .lt. -10.) goto 20
      call splint(eta,f12,f12a,n,y,f1,klo,khi)
      GO TO 100
 10   X2=y**(-2)
      F1=A(6)*y*SQRT(y)*th*(5+(A(1)+(3*A(2)+7*X2*A(3))*X2)*x2)
      GO TO 100
 20   F0=DEXP(y)
      F1=A(7)*th*F0*(2-(4*A(4)-(6*A(5)-.25*F0)*F0)*F0)
 100  F_1_2=F1
 999  RETURN
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       F_3_2
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         11/29/89
C
C    CALL LINE:    F_3_2(Y)      (3/2th Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       3/2th Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION F_3_2(y)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(n=201)
      DIMENSION A(7),eta(n),f32(n),f32a(n),f12(n),fr(n),f12a(n),fra(n),
     > fia(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
      DATA A,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0
     1,6.41500299D-02,.4D0,1.32934039D0,1,n/
      IF(y .gt. 30.) goto 10
      if(y .lt. -10.) goto 20
      call splint(eta,f32,f32a,n,y,f1,klo,khi)
      GO TO 100
 10   X2=y**(-2)
      F1=A(6)*SQRT(y)*(1./X2+A(1)-(A(2)+X2*A(3))*X2)
      GO TO 100
 20   F0=DEXP(y)
      F1=A(7)*F0*(1.-(A(4)-(A(5)-.03125*F0)*F0)*F0)
 100  F_3_2=F1
 999  RETURN
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       FHALF
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         3/10/90
C
C    CALL LINE:    FHALF(Y)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       Ratio of 1/2th Fermi Integral to the -1/2th Fermi
C                  Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION FHALF(y)
      IMPLICIT real*8(a-h,o-z)
      parameter(n=201)
      DIMENSION A(7),f12(n),fia(n),eta(n),f32(n),fr(n),f32a(n),
     > fra(n),f12a(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
      DATA A,th,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0
     1,6.41500299D-02,.4D0,1.32934039D0,.3333333333333d0,1,n/
      IF(y .gt. 30.) goto 10
      if(y .lt. -10.) goto 20
      call splint(eta,fr,fra,n,y,f1,klo,khi)
      GO TO 100
 10   X2=y**(-2)
      F1=y*th*(1.+(.2*A(1)+(.6*A(2)+1.4*X2*A(3))*X2)*x2)
     > /(1.-(.2*th*a(1)-(a(2)-4.2*x2*a(3))*x2)*x2)
      GO TO 100
 20   F0=EXP(y)
      F1=(1.-(2*a(4)-(3*a(5)-.125*f0)*f0)*f0)/
     > (2.-(8*a(4)-(18*a(5)-f0)*f0)*f0)
 100  FHALF=F1
 999  RETURN
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       FINV12
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         11/29/89
C
C    CALL LINE:    FINV12(Y)      (Inverse of the 1/2th Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       Inverse of Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION FINV12(y)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(n=201)
      DIMENSION AI(8),f12(n),fia(n),eta(n),f32(n),fr(n),f32a(n),
     > fra(n),f12a(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
      DATA AI,klo,khi/-.822467032D0,-1.21761363D0,-9.16138616D0,
     1.398942281D0,.0732748216D0,-1.310707D0,1.12837917D0,
     28.2810645D-3,1,n/
      if(y .gt. 109.695) goto 10
      if(y .lt. 4.0234e-5) goto 20
      call splint(f12,eta,fia,n,y,f1,klo,khi)
      GO TO 100
 10   X2=(1.5*y)**(.666666667)
      X4=1./(X2*X2)
      F1=X2*(1.+(AI(1)+(AI(2)+AI(3)*X4)*X4)*X4)
      GO TO 100
 20   F1=LOG(AI(7)*MAX(y,1.D-20)*(1.+(AI(4)+(AI(5)+AI(8)*y)*y)*y))
 100  finv12=F1
 999  RETURN
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       GETFNM.FOR
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         8/5/89
C
C    PURPOSE:      OBTAINS A FILE NAME FROM USER
C
C    CALL LINE:    CALL GETFNM(FNAME,FNML,PROMPT,PROMTL)
C
C    INPUTS:       PROMPT = STRING TO POMPT USER WITH (C*60)
C                  PROMTL = LENGTH OF PROMPT STRING (I)
C
C    OUTPUTS:      FNAME = FILE NAME (C*60)
C                  FNML = FILE NAME LENGTH (I)
C*************************************************************************
C
      SUBROUTINE GETFNM(FNAME,FNML,PROMPT,PROMTL)
C
      IMPLICIT NONE
C
      INTEGER FNML, PROMTL
      CHARACTER*60 FNAME, PROMPT
C
C                       Local variables
C
      INTEGER STDOUT, STDIN, I
      DATA STDOUT/6/, STDIN/5/
C
C                       Prompt user for file name
C
      WRITE(STDOUT,'(T2,A,$)') PROMPT(1:PROMTL)
      READ(STDIN,'(A)') FNAME
C
C                        Figure out input file name length
      DO 10 I=1,20,1
        IF(FNAME(I:I).EQ.' ') GOTO 20
 10   CONTINUE
C
 20   CONTINUE
      FNML = I-1
C
C
 999  RETURN
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         INVEOS
C    MODULE:       INVEOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         5/23/90
C                  Bug fixed on (5/24/90) (affected only performance
C                  of code NOT the results!)
C
C
C    CALL LINE:    CALL INVEOS(INPVAR,T_OLD,YE,BRYDNS,IFLAG,EOSFLG,XPREV)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  T_OLD = INITIAL GUESS AT THE TEMPERATURE
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C                  IFLAG = 1 --> INPVAR IS TEMPERATURE
C                          2 --> INPVAR IS INTERNAL ENERGY
C                          3 --> INPVAR IS ENTROPY (NOT IMPLEM)
C
C    OUTPUTS       EOSFLG = 1 --> "NO NUCLEI" EOS
C                           2 --> GENERAL EOS
C                           3 --> BULK EOS FOR DENSITIES ABOVE NUCLEAR
C                  XPREV = UNUSED
C                  P_PREV = PREVIOUS VALUE OF PROTON DENSITY
C
C
C
C
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE INVEOS(INPVARP,T_OLD,YEP,BRYDNSP,IFLAGP,EOSFLGP,
     1                  FORFLGP,SF,XPREVP,P_PREVP)
C
C
C
 
C
      USE eos_m4c_module
C
      IMPLICIT NONE
C
C                         Local variables
C
      DOUBLE PRECISION INP_V, INP_VO, INP_VN, UFTN, DUFTN, DT
      DOUBLE PRECISION T_OLD, T_NEW, T_TEMP, T_LB, T_UB, PERDIF
      INTEGER LOOP, SF, NEW_F
C
      INTEGER IFLAGP,EOSFLGP,FORFLGP
      DOUBLE PRECISION INPVARP(4),YEP,BRYDNSP,XPREVP,P_PREVP
C
C
      INPVAR(1)      = INPVARP(1)
      INPVAR(2)      = INPVARP(2)
      INPVAR(3)      = INPVARP(3)
      INPVAR(4)      = INPVARP(4)
      YE             = YEP
      BRYDNS         = BRYDNSP
      IFLAG          = IFLAGP
      EOSFLG         = EOSFLGP
      FORFLG         = FORFLGP
      XPREV          = XPREVP
      P_PREV         = P_PREVP
C
      RSFLAG = 1
C                         Input is the temperature; call the EOS
C                         normally and then return
      IF(IFLAG.EQ.1) THEN
        CALL EOS_M4C(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,
     1 XPREV,P_PREV)
       T_OLD = INPVAR(1)
        INPVARP(1)   = INPVAR(1)
        INPVARP(2)   = INPVAR(2)
        INPVARP(3)   = INPVAR(3)
        INPVARP(4)   = INPVAR(4)
        YEP          = YE
        BRYDNSP      = BRYDNS
        IFLAGP       = IFLAG
        EOSFLGP      = EOSFLG
        FORFLGP      = FORFLG
        XPREVP       = XPREV
        P_PREVP      = P_PREV
        RETURN
      ENDIF
C
C
C                         The input variable must be the internal
C                         energy so calc the internal energy for
C                         the initial guess at the temperature
        INP_V = INPVAR(1)
C
        T_LB = 0.05
        T_UB = 50.0
C
        INPVAR(1) = T_OLD
        CALL EOS_M4C(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,
     1 XPREV,P_PREV)
CCC      CALL EOS_M1D(T_OLD,YE,BRYDNS,1,EOSFLG,XPREV,P_PREV)
C
C                         Save the value of the internal energy
      IF(IFLAG.EQ.2) THEN
        INP_VO = UTOT
      ELSEIF(IFLAG.EQ.3) THEN
        INP_VO = STOT
      ENDIF
C
C
C                         Tweak the initial guess slightly so as to
C                         get a new value of the internal energy
C
      T_NEW = 1.1*T_OLD
C
      NEW_F = 1
C
      DO 20 LOOP=1,50,1
C
        INPVAR(1) = T_NEW
        CALL EOS_M4C(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,
     1 XPREV,P_PREV)
CCC        CALL EOS_M1D(T_NEW,YE,BRYDNS,1,EOSFLG,XPREV,P_PREV)
C
C
        IF(SF.NE.1.AND.NEW_F.EQ.1) THEN
cc          WRITE(*,*) 'INVEOS: EOS fatally failed at try:'
cc          WRITE(*,*) T_NEW,BRYDNS,YE
          T_NEW = T_NEW-0.25
          T_LB = DMIN1(T_LB,T_NEW-1.0D-1)
          GOTO 20
        ELSEIF(SF.NE.1) THEN
           DT = 0.5*DT
           T_NEW = T_NEW+DT
        ELSE
C
          NEW_F = 0
C
C                         Save this value of the internal energy too
          IF(IFLAG.EQ.2) THEN
            INP_VN = UTOT
          ELSEIF(IFLAG.EQ.3) THEN
            INP_VN = STOT
          ENDIF
C
          IF(INP_VN.LT.INP_V) THEN
            T_LB = T_NEW
c            write(*,*) 'l @ ',t_new,inp_vn,inp_v
          ELSEIF(INP_VN.GT.INP_V) THEN
c            write(*,*) 'u @ ',t_new,inp_vn,inp_v
            T_UB = T_NEW
          ENDIF
        ENDIF
C
        UFTN = INP_VN-INP_V
C
        IF(LOOP.LT.20) THEN
C                         This is the function to be zeroed by the
C                         Newton-Raphson iteration
C
C                         Numerical derivative of the above function
C                         w.r.t. the temperature
CC          DUFTN = ((INP_VN-INP_VO)/(T_NEW-T_OLD))+1.0D-15
C
C                         Analytic derivatives
          IF(IFLAG.EQ.2) THEN
            DUFTN = DUDT
          ELSEIF(IFLAG.EQ.3) THEN
            DUFTN = DSDT
          ENDIF
C
C                         Estimated correction to temperature
          DT = UFTN/DUFTN
C
C
 10       CONTINUE
C                         Temporarily store the new temperature
          T_TEMP = T_NEW-DT
C
C                         Is the new temp within a valid range?
          IF((T_TEMP.GT.T_LB).AND.(T_TEMP.LT.T_UB)) THEN
            T_OLD = T_NEW
            INP_VO = INP_VN
            T_NEW = T_TEMP
          ELSE
C                         If not cut the step size in half & try again
            DT = 0.5*DT
            IF(T_TEMP.EQ.T_NEW) THEN
              T_OLD = T_NEW
              T_NEW = 0.5*(T_LB+T_UB)
              DT = T_NEW-T_OLD
            ELSE
              GOTO 10
            ENDIF
          ENDIF
C
        ELSE
C
          T_OLD = T_NEW
          T_NEW = 0.5*(T_LB+T_UB)
          DT = T_NEW-T_OLD
C
        ENDIF
C
C                         If relative change in T is less than 1.0e-5
C                         then quit out of loop
        IF(ABS(DT/T_NEW).LT.1.0D-5) GOTO 30
C
C                End of the Do loop
 20   CONTINUE
C
C                Didn't meet convergence criterion
      SF = 0
      WRITE(*,*) ' INVERSION OF EOS FAILED TO CONVERGE',DT,T_NEW
C
C                Met the convergence criterion!!!
 30   CONTINUE
C
c                This stuff is commented out for speed reasons; it
c                virtually never gets tripped anyway
cc      INPVAR(1) = T_OLD
cc      CALL EOS_M4C(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,XPREV,P_PREV)
cc      IF(IFLAG.EQ.2) THEN
cc        PERDIF = INP_V-UTOT
cc      ELSE
cc        PERDIF = INP_V-STOT
cc      ENDIF
cc      IF(ABS(PERDIF).GT.1.0D-4) THEN
cc        WRITE(*,*) 'INVEOS: FAILURE',INP_V,STOT
cc        write(*,*) uftn,dt,loop
cc        WRITE(*,*) 'TRY:',T_NEW,BRYDNS,YE
cc        SF = 0
cc        RETURN
cc      ENDIF
C
C                Return this value for T
      INPVAR(1) = INP_V
      T_OLD = T_NEW
C
C                Time to call it quits!
  999 CONTINUE
C
      INPVARP(1)     = INPVAR(1)
      INPVARP(2)     = INPVAR(2)
      INPVARP(3)     = INPVAR(3)
      INPVARP(4)     = INPVAR(4)
      YEP            = YE
      BRYDNSP        = BRYDNS
      IFLAGP         = IFLAG
      EOSFLGP        = EOSFLG
      FORFLGP        = FORFLG
      XPREVP         = XPREV
      P_PREVP        = P_PREV
C
      RETURN
C
      END
C
C

!///////////////////////////////////////////////////////////////////////

c----------------------------------------------------------------------c
c                                                                      c
c    File:         lubksb                                              c
c    Module:       lubksb                                              c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         8/22/00                                             c
c                                                                      c
c    Purpose:                                                          c
c      Solves the set of of n linear equations ax = b. Here a is       c
c       input, not as the matrix a but rather as its LU decomposition, c
c       determined by subroutine ludcmp. indx is input as the          c
c       permutation vector returned by ludcmp. b is input as the       c
c       right-hand side vector, and returns with the solution vector   c
c       x. a, n, np, and indx are not modified by this rotine and      c
c       can be left in place for successive calls with different       c
c       right-hand sides b. This routine takes into account the        c
c       possibility that b will begin with many zero elements, so it   c
c       is efficient for use in matrix inversion.                      c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c   a         : nxn LU decomposed matrix of coefficients               c
c   n         : dimension of a                                         c
c   np        : physical dimension of a                                c
c   indx      : record of row permutation                              c
c   b         : constant vector                                        c
c                                                                      c
c    Output arguments:                                                 c
c   b         : solution vector                                        c
c                                                                      c
c    Variables that must be passed through common:                     c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      none                                                            c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine lubksb(a,n,np,indx,b)
c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer i,ii,indx,j,ll,n,np
      double precision a,b,sum,zero
c                                                                      c
c                                                                      c
      dimension a(np,np),indx(n),b(n)
c                                                                      c
c                                                                      c
      data zero /0.0d+00/
c----------------------------------------------------------------------c
c        When ii is set to a positive value, it will become the index  c
c         of the first nonvanishing element of b. We now do the        c
c         forward substitution. The only new wrinkle is to unscramble  c
c         the permutation as we go.                                    c
c----------------------------------------------------------------------c
      ii                = 0
      do 1000 i = 1,n
       ll               = indx(i)
       sum              = b(ll)
       b(ll)            = b(i)
       if ( ii .ne. 0 ) then
        do 100 j = ii,i-1
         sum            = sum - a(i,j) * b(j)
  100   continue
       else if ( sum .ne. zero ) then
        ii              = i
       end if ! ii ne 0
       b(i)             = sum
 1000 continue
c----------------------------------------------------------------------c
c        Now we do the back substitution.                              c
c----------------------------------------------------------------------c
      do 2000 i = n,1,-1
       sum              = b(i)
       if ( i .lt. n ) then
       do 200 j = i+1,n
        sum             = sum - a(i,j) * b(j)
  200  continue
       end if
       b(i)             = sum/a(i,i)
 2000 continue
c                                                                      c
c                                                                      c
      return
      end
      
!///////////////////////////////////////////////////////////////////////

c----------------------------------------------------------------------c
c                                                                      c
c    File:         ludcmp                                              c
c    Module:       ludcmp                                              c
c    Type:         Subprogram                                          c
c    Author:       S. W. Bruenn, Dept of Physics, FAU,                 c
c                  Boca Raton, FL 33431-0991                           c
c                                                                      c
c    Date:         8/21/00                                             c
c                                                                      c
c    Purpose:                                                          c
c      Given an nxn matrix a, with physical dimension np, this         c
c       replaces it by the LU decomposition of a rowwise permutation   c
c       of itself. A is output arranged as                             c
c                                                                      c
c              b11  b12  b13  b14                                      c
c              a21  b22  b23  b24                                      c
c              a31  a32  b33  b34                                      c
c              a41  a42  a43  b44                                      c
c                                                                      c
c       indx is an output vector which records the row permutation     c
c       effected by the partial pivoting; d is output + or - 1         c
c       depending on whether the number of row interchanges was        c
c       even or odd, respectively. This routine is used in combina-    c
c       tion with lubksb to solve linear equations or invert a         c
c       matrix.                                                        c
c                                                                      c
c    Subprograms called:                                               c
c      none                                                            c
c                                                                      c
c    Input arguments:                                                  c
c   a         : nxn matrix                                             c
c   n         : dimension of a                                         c
c   np        : physical dimension of a                                c
c                                                                      c
c    Output arguments:                                                 c
c   a         : LU decomposition of a                                  c
c   indx      : record of row permutation                              c
c   d         : odd or even permutation number index                   c
c                                                                      c
c    Variables that must be passed through common:                     c
c      none                                                            c
c                                                                      c
c    Include files:                                                    c
c      none                                                            c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine ludcmp(a,n,np,indx,d)

c                                                                      c
c                                                                      c
      implicit none
      save
c                                                                      c
c        Local variables (not in common blocks)                        c
      integer i,imax,indx,j,j1,j2,k,n,np,nmax
      double precision a,aamax,d,dum,one,sum,tiny,vv,zero
c                                                                      c
c                                                                      c
      parameter (nmax = 100,tiny = 1.d-20)
c                                                                      c
c                                                                      c
      dimension a(np,np),indx(n),vv(nmax)
c                                                                      c
c                                                                      c
      data zero /0.0d+00/
      data one  /1.0d+00/

 1001 FORMAT (81(1pe11.3))

c----------------------------------------------------------------------c
c        Loop over rows to get the implicit scaling information.       c
c----------------------------------------------------------------------c
      d             = one
      do 1000 i = 1,n
       aamax        = zero
       do 100 j = 1,n
        if ( dabs(a(i,j)) .gt. aamax ) aamax = dabs(a(i,j))
  100  continue
       IF ( aamax .eq. zero ) THEN
         WRITE (6,*) 'Singular matrix.'
         DO j1 = 1,n
           WRITE (6,1001) (a(j1,j2),j2 = 1,n)
         END DO
         STOP
       END IF
       vv(i)        = one/aamax
 1000 continue
c----------------------------------------------------------------------c
c        Loop over columns of Crout's method.                          c
c----------------------------------------------------------------------c
      do 2000 j = 1,n
       if ( j .gt. 1 ) then
        do 200 i = 1,j-1
         sum        = a(i,j)
         if ( i .gt. 1 ) then
          do 20 k = 1,i-1
           sum      = sum - a(i,k) * a(k,j)
   20     continue
          a(i,j)    = sum
         end if ! i > 1
  200   continue
       end if ! j > 1
c----------------------------------------------------------------------c
c        Initialize the search for the largest pivot.                  c
c----------------------------------------------------------------------c
       aamax        = zero
       do 210 i = j,n
        sum         = a(i,j)
        if ( j .gt. 1 ) then
         do 21 k = 1,j-1
          sum       = sum - a(i,k) * a(k,j)
   21    continue
         a(i,j)     = sum
        end if ! j > 1
        dum         = vv(i) * dabs(sum)
        if ( dum .ge. aamax ) then
         imax       = i
         aamax      = dum
        end if ! dum ge aamax
  210  continue
c----------------------------------------------------------------------c
c        Do we need to interchange rows?                               c
c----------------------------------------------------------------------c
       if ( j .ne. imax ) then
        do 220 k = 1,n
         dum        = a(imax,k)
         a(imax,k)  = a(j,k)
         a(j,k)     = dum
  220   continue
        d           = -d
        vv(imax)    = vv(j)
       end if ! j ne imax
       indx(j)      = imax
c----------------------------------------------------------------------c
c        If the pivot element is zero the matrix is singular (at       c
c         least to the precision of the algorithm). For some applica-  c
c         tions on singular matrices, it is desirable to substitute    c
c         tiny for zero.                                               c
c----------------------------------------------------------------------c
       if ( a(j,j) .eq. zero ) a(j,j) = tiny
c----------------------------------------------------------------------c
c        Now, finally, divide by the pivot element.                    c
c----------------------------------------------------------------------c
       if ( j .ne. n ) then
        dum         = one/a(j,j)
        do 230 i = j+1,n
         a(i,j)     = a(i,j) * dum
  230   continue
       end if ! j ne n
c----------------------------------------------------------------------c
c        Go back for the next column in the reduction.                 c
c----------------------------------------------------------------------c
 2000 continue
c                                                                      c
c                                                                      c
      return
      end

!///////////////////////////////////////////////////////////////////////

C***********************************************************************
C
C    MODULE:       MATADD
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Adds two N by M Matrices
C
C    CALL LINE:    CALL MATINV(A,B,C,N,M)
C
C    INPUTS:       A,B= Arrays to be added  (D)
C                  N = Number of rows in arrays (I)
C                  M = Number of columns in arrays (I)
C
C    OUTPUTS:      C = Array containing A+B (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATADD(A,B,C,N,M)
C
      IMPLICIT NONE
C
      INTEGER N,M
      DOUBLE PRECISION A(N,M), B(N,M), C(N,M)
C
C
C                 Local variables
C
      INTEGER I, J
C
      DO 20 J=1,M,1
        DO 10 I=1,N,1
          C(I,J) = A(I,J)+B(I,J)
 10     CONTINUE
 20   CONTINUE
C
 999  RETURN
C
      END

!///////////////////////////////////////////////////////////////////////

C
C***********************************************************************
C
C    MODULE:       MATCOP
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Copy one N by M Matrix into another
C
C    CALL LINE:    CALL MATCOP(A,B,N,M)
C
C    INPUTS:       A = Array to be copied  (D)
C                  N = Number of rows in arrays (I)
C                  M = Number of columns in arrays (I)
C
C    OUTPUTS:      C = Array to be copied into (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATCOP(A,B,N,M)
C
      IMPLICIT NONE
C
      INTEGER N,M
      DOUBLE PRECISION A(N,M), B(N,M)
C
C
C                 Local variables
C
      INTEGER I, J
C
      DO 20 J=1,M,1
        DO 10 I=1,N,1
          B(I,J) = A(I,J)
 10     CONTINUE
 20   CONTINUE
C
 999  RETURN
C
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         DMATRIX.FOR
C
C***********************************************************************
C
C    MODULE:       MATINV
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Inverts a N by N Matrix
C
C    CALL LINE:    CALL MATINV(A,AINV,N)
C
C    INPUTS:       A = Array to be inverted  (D)
C                  N = dimesion of arrays (I)
C
C    OUTPUTS:      AINV = Inverse of A (D)
C
C    CALLS :       Numerical recipes routines LUDCMP, LUBKSB
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATINV(A,AINV,N)
C
      IMPLICIT NONE
C
      INTEGER N
      DOUBLE PRECISION A(N,N), AINV(N,N)
C
C
C                 Local variables
C
      INTEGER NPHYS, I, J
      PARAMETER(NPHYS=10)
      DOUBLE PRECISION TEMP(NPHYS,NPHYS), Y(NPHYS,NPHYS), D
      INTEGER INDEX(NPHYS)
C
C                 Make a copy of the array, and initialize
C                 the indentity matrix
      DO 20 J=1,N,1
        DO 10 I=1,N,1
          Y(I,J) = 0.0
          TEMP(I,J) = A(I,J)
 10     CONTINUE
        Y(J,J) = 1.0
 20   CONTINUE
C
C
C                 LU decompose the matrix
      CALL LUDCMP(TEMP,N,NPHYS,INDEX,D)
C
C
C                 Back substitute to get inverse
      DO 30 J=1,N,1
        CALL LUBKSB(TEMP,N,NPHYS,INDEX,Y(1,J))
 30   CONTINUE
C
C
C                 Copy temporary array into the inverse array
      DO 50 J=1,N,1
        DO 40 I=1,N,1
          AINV(I,J) = Y(I,J)
 40     CONTINUE
 50   CONTINUE
C
C
 999  RETURN
      END

!///////////////////////////////////////////////////////////////////////

C
C***********************************************************************
C
C    MODULE:       MATMUL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Multiplies two Matrices (LxM)x(MxN)
C
C    CALL LINE:    CALL MATMUL(A,B,C,L,M,N)
C
C    INPUTS:       A,B= Arrays to be added  (D)
C                  L,M,N = Dimensions of arrays (I)
C
C    OUTPUTS:      C = Array containing A x B (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATMUL(A,B,C,L,M,N)
C
      IMPLICIT NONE
C
      INTEGER L, M, N
      DOUBLE PRECISION A(L,M), B(M,N), C(L,N)
C
C
C                 Local variables
C
      INTEGER I, J, K
      DOUBLE PRECISION SUM
C
C                 Loop over all elements of the array
      DO 30 I=1,L,1
        DO 20 J=1,N,1
C
C                 Initialize SUM for a new element
          SUM = 0.0
C                 Calculate (i,j)th element
          DO 10 K=1,M,1
            SUM = SUM+A(I,K)*B(K,J)
 10       CONTINUE
          C(I,J) = SUM
C
 20     CONTINUE
 30   CONTINUE
C
 999  RETURN
C
      END

!///////////////////////////////////////////////////////////////////////

C
C***********************************************************************
C
C    MODULE:       MATSCL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Multiply a N by M Matrix by a scalar
C
C    CALL LINE:    CALL MATSCL(A,SCALAR,B,N,M)
C
C    INPUTS:       A = Array to be scaled  (D)
C                  SCALAR = Constant to multiply matrix by (D)
C                  N = Number of rows in arrays (I)
C                  M = Number of columns in arrays (I)
C
C    OUTPUTS:      B = Array containing SCALAR x A (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATSCL(A,SCALAR,B,N,M)
C
      IMPLICIT NONE
C
      INTEGER N, M
      DOUBLE PRECISION A(N,M), B(N,M), SCALAR
C
C
C                 Local variables
C
      INTEGER I, J
C
C                 Loop over all elements of the array
      DO 20 J=1,M,1
        DO 10 I=1,N,1
C
          B(I,J) = SCALAR*A(I,J)
C
 10     CONTINUE
 20   CONTINUE
C
 999  RETURN
C
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         MAXWEL.FOR
C
C***********************************************************************
C
C    MODULE:       MAXWEL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         7/13/90
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU
C
C
C    CALL LINE:    CALL MAXWEL(INPVAR,YE,BRYDNS)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:
C
C
C
C
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE MAXWEL(INPVARP,YEP,BRYDNSP,XPREVP,P_PREVP,SSFLAGP)
C
C
      USE el_eos_module
      USE eos_m4c_module
      USE maxwel_module
C
      IMPLICIT NONE
C
      DOUBLE PRECISION OUTVAR(4)
      DOUBLE PRECISION DPN_DT, DPN_DN, DPN_DY
      DOUBLE PRECISION DSN_DT, DSN_DN, DSN_DY
      DOUBLE PRECISION DSB_DT, DSB_DN, DSB_DY
      DOUBLE PRECISION DMU_DT, DMU_DN, DMU_DY
      DOUBLE PRECISION DPHADT, DPHADY, DELDNS
      DOUBLE PRECISION N_XH, N_XA, N_XN, N_XP, B_XA, B_XN, B_XP
C
      INTEGER SSFLAGP
      DOUBLE PRECISION INPVARP(4),YEP,BRYDNSP,XPREVP,P_PREVP,BRYDNS_T
C
C
      INPVAR(1)      = INPVARP(1)
      INPVAR(2)      = INPVARP(2)
      INPVAR(3)      = INPVARP(3)
      INPVAR(4)      = INPVARP(4)
      YE             = YEP
      BRYDNS         = BRYDNSP
      XPREV          = XPREVP
      P_PREV         = P_PREVP
      SSFLAG         = SSFLAGP
C
C
C                   Set the temperature
      T = INPVAR(1)
C
C
C                   Calculate and save chem. pot. and thermodynamic
C                   quantaties from low end of two phase region
      BRYDNS_T = BRYDNS
      CALL NUCEOS(INPVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG)
      BRYDNS   = BRYDNS_T 
C
C
C
C                    If the nuclear EOS failed and the reset flag is set
C                    then reset the initial guesses and try again
      IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN
        BRYDNS_T = BRYDNS
        CALL RESET(INPVAR,YE,LOWDNS,OUTVAR)
        BRYDNS   = BRYDNS_T 
        OUTVAR(1) = INPVAR(1)
        BRYDNS_T = BRYDNS
        CALL NUCEOS(OUTVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG)
        BRYDNS   = BRYDNS_T 
C
C
C                    Make a last ditch effort at convergence
        IF(SSFLAG.NE.1) THEN
          OUTVAR(2) = 0.155
          OUTVAR(3) = -15.0
          OUTVAR(4) = -20.0
          BRYDNS_T = BRYDNS
          CALL NUCEOS(OUTVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG)
          BRYDNS   = BRYDNS_T 
        ELSE
          INPVAR(2) = OUTVAR(2)
          INPVAR(3) = OUTVAR(3)
          INPVAR(4) = OUTVAR(4)
        ENDIF
C
      ENDIF
C
C
C
C
      PRLOW = PTOT-PPRESS
      S_LOW = STOT-PS
      F_LOW = FTOT-PF
      MUTLOW = (1.0-YE)*MUN+YE*(MUPROT+MUSUBE)
      MUELOW = MUSUBE
      MUHLOW = MUHAT
C
      DPN_DT = DPDT
      DPN_DN = DPDN
      DPN_DY = DPDY
C
      DMU_DT = DMUDT
      DMU_DN = DMUDN
      DMU_DY = DMUDY
C
      DSN_DT = DSDT-DPSDT
      DSN_DN = DSDN
      DSN_DY = DSDY
C
      N_XH = XH
      N_XA = XALFA
      N_XP = XPROT
      N_XN = XNUT
C
C
      IF(SSFLAG.NE.1) THEN
        WRITE(*,*) 'MAXWEL:  Nuclear EOS failed at try:'
        WRITE(*,*) T,LOWDNS,YE
        WRITE(*,*) INPVAR
        GOTO 999
      ENDIF
C                   Calculate and save chem. pot. and thermodynamic
C                   quantaties from high end of two phase region
      BRYDNS_T = BRYDNS
      CALL ALFEOS(INPVAR,YE,HIDNS,P_PREV,SSFLAG)
      BRYDNS   = BRYDNS_T 
C
      PRHI = PTOT-PPRESS
      S_HI = STOT-PS
      F_HI = FTOT-PF
      MUTHI = (1.0-YE)*MUN+YE*(MUPROT+MUSUBE)
      MUEHI = MUSUBE
      MUHHI = MUHAT
C
C
      DSB_DT = DSDT-DPSDT
      DSB_DN = DSDN
      DSB_DY = DSDY
C
C
      B_XA = XALFA
      B_XP = XPROT
      B_XN = XNUT
C
C
      IF(SSFLAG.NE.1) THEN
        WRITE(*,*) 'MAXWEL:  Alfa EOS failed at try:'
        WRITE(*,*) T,HIDNS,YE
        WRITE(*,*) INPVAR
        GOTO 999
      ENDIF
C
C                   Calculate "average" chem. pot. and pressure
C                   in order to avoid numerical problems
      MUTILD = (MUTLOW+MUTHI)/2.0
      PRTILD = (PRLOW+PRHI)/2.0
C
C                   Calculate phase fraction
      PHASEF = (BRYDNS-LOWDNS)/(HIDNS-LOWDNS)
C
C
C                   Electron number density
      NSUBE = BRYDNS*YE
C
C                   Call electron EOS to determine the
C                   electron chemical potential
      CALL EL_EOS(T,YE,BRYDNS)
C
C
      MUHAT = MUSUBE+(1.0-PHASEF)*(MUHLOW-MUELOW)+PHASEF*(MUHHI-MUEHI)
C
      MUN = MUTILD+YE*(MUHAT-MUSUBE)
C
      MUPROT = MUN-MUHAT
C
C                   Calculate thermodynamic quantities
C
      STOT = ((1.0-PHASEF)*S_LOW*LOWDNS+PHASEF*S_HI*HIDNS)/BRYDNS+PS
C
      FTOT = (LOWDNS*F_LOW+MUTILD*(BRYDNS-LOWDNS))/BRYDNS+PF
C
      UTOT = FTOT+T*STOT+PU
C
      PTOT = PRTILD+PPRESS
C
C
      XH = (1.0-PHASEF)*N_XH
      XALFA = (1.0-PHASEF)*N_XA
      XNUT = (1.0-PHASEF)*N_XN
      XPROT = (1.0-PHASEF)*N_XP
      XALFA2 = PHASEF*B_XA
      XNUT2 = PHASEF*B_XN
      XPROT2 = PHASEF*B_XP
C
C
C
C
      DELDNS = HIDNS-LOWDNS
C
C
      DPHADT = ((BRYDNS-LOWDNS)/DELDNS**2-1.0/DELDNS)*DNL_DT-
     1    ((BRYDNS-LOWDNS)/DELDNS**2)*DNH_DT
C
      DPDT = DPN_DT+DPN_DN*DNL_DT
      DMUDT = DMU_DT+DMU_DN*DNL_DT
      DSDT = (1.0-PHASEF)*LOWDNS*(DSN_DT+DSN_DN*DNL_DT)/BRYDNS+
     2 (1.0-PHASEF)*S_LOW*DNL_DT/BRYDNS-LOWDNS*S_LOW*DPHADT/BRYDNS+
     3    (DPHADT*S_HI*HIDNS+PHASEF*DNH_DT*S_HI+
     4    PHASEF*HIDNS*(DSB_DT+DSB_DN*DNH_DT))/BRYDNS+DPSDT
      DUDT = DMUDT-DPDT/BRYDNS+STOT+T*DSDT
C
C
      DPDN = 0.0
      DMUDN = 0.0
      DSDN = -DPDT/BRYDNS**2
      DUDN = (LOWDNS*(MUTILD-FTOT)/BRYDNS**2)+T*DSDN
C
C
      DPHADY = ((BRYDNS-LOWDNS)/DELDNS**2-1.0/DELDNS)*DNL_DY-
     1    ((BRYDNS-LOWDNS)/DELDNS**2)*DNH_DY
C
      DPDY = DPN_DY+DPN_DN*DNL_DY
      DMUDY = DMU_DY+DMU_DN*DNL_DY
      DSDY = (1.0-PHASEF)*LOWDNS*(DSN_DY+DSN_DN*DNL_DY)/BRYDNS+
     2 (1.0-PHASEF)*S_LOW*DNL_DY/BRYDNS-LOWDNS*S_LOW*DPHADY/BRYDNS+
     3    (DPHADY*S_HI*HIDNS+PHASEF*DNH_DY*S_HI+
     4    PHASEF*HIDNS*(DSB_DY+DSB_DN*DNH_DY))/BRYDNS
      DUDY = DMUDY-DPDY/BRYDNS+T*DSDY
C
C
C
C
C             Adiabatic index
C             (Note that the first term vanishes in this expression)
      GAM_S = T*(DPDT**2)/(BRYDNS*PTOT*DUDT)
C
  999 CONTINUE
C
C
      INPVARP(1)     = INPVAR(1)
      INPVARP(2)     = INPVAR(2)
      INPVARP(3)     = INPVAR(3)
      INPVARP(4)     = INPVAR(4)
      YEP            = YE
      BRYDNSP        = BRYDNS
      XPREVP         = XPREV
      P_PREVP        = P_PREV
      SSFLAGP        = SSFLAG
C
C
      RETURN
C
C
      END

!///////////////////////////////////////////////////////////////////////

C***********************************************************************
C
C    MODULE:       MV_MUL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/3/90
C
C    PURPOSE:      Multiplies a Matrix times a vector (NxN)x(N)
C
C    CALL LINE:    CALL MV_MUL(A,V,RV,N)
C
C    INPUTS:       A = Array to be multiplied  (D)
C                  V = Vector to be multiplied (D)
C                  N = Dimensions of arrays & vector (I)
C
C    OUTPUTS:      RV = resultant vector (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MV_MUL(A,V,RV,N)
C
      IMPLICIT NONE
C
      INTEGER N
      DOUBLE PRECISION A(N,N), V(N), RV(N)
C
C
C                 Local variables
C
      INTEGER I, J
      DOUBLE PRECISION SUM
C
C                 Loop over all elements of the array
      DO 20 I=1,N,1
C
C                 Initialize SUM for a new element
        SUM = 0.0
C                 Calculate (i)th element
        DO 10 J=1,N,1
          SUM = SUM+A(I,J)*V(J)
 10     CONTINUE
        RV(I) = SUM
C
 20   CONTINUE
C
 999  RETURN
C
      END

!///////////////////////////////////////////////////////////////////////

Cnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnu
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         NUCEOS.FOR
C
C***********************************************************************
C
C    MODULE:       NUCEOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         7/13/90 Modified from model 1-d
C
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C    CALL LINE:    CALL NUCEOS(INPVAR,YE,BRYDNS,X_PREV,SSFLAG)
C
C    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C    OUTPUTS:      XPREV = PREVIOUS VALUE OF X (MUST BE SUPPLIED ON
C                  FIRST CALL)
C                  SSFLAG = SUCCESS FLAG 0 --> FAILURE
C                                        1 --> SUCCESS
C
C
C
C    INCLUDE FILES:  EOS_M4C.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE NUCEOS(INPVARP,YEP,BRYDNSP,XPREVP,P_PREVP,SSFLAGP)
C
C
      USE el_eos_module
      USE eos_m4c_module
C
      IMPLICIT NONE
C
C
C                       Function type declarations
C
      DOUBLE PRECISION F_1_2, F_3_2, FINV12, FHALFI, FHALFO
      double precision fhalf
C
      DOUBLE PRECISION ZNG, ZPG
      INTEGER TCFLAG, ftflag
C
      INTEGER KKI,LLI
      DOUBLE PRECISION RESULT(5), R_CHECK(5)
      double precision a_tmp(5,5)
      DOUBLE PRECISION NI_MIN
C
      integer cflag, schflg
      double precision dtst1, dtst2
      double precision break, dnsi, dtmp8
      double precision dtmp1,dtmp2,dtmp3,dtmp4,dtmp5,dtmp6,dtmp7
C
      INTEGER SSFLAGP
      DOUBLE PRECISION INPVARP(4),YEP,BRYDNSP,XPREVP,P_PREVP
cc      double precision tbsph, tbph, tbnh, tbspl, tbpl, tbnl
cc      double precision dbspdx, dbpdx, dbndx, dbspdu, dbpdu, dbndu
cc      double precision tsgl, tsgh, thl, thh, dsgdx, dhfdx, ds2dx,dzdx
cc      double precision dpt1dx, dpt2dx
c
      INCLUDE 'force.inc'
C
C
      INPVAR(1)      = INPVARP(1)
      INPVAR(2)      = INPVARP(2)
      INPVAR(3)      = INPVARP(3)
      INPVAR(4)      = INPVARP(4)
      YE             = YEP
      BRYDNS         = BRYDNSP
      XPREV          = XPREVP
      P_PREV         = P_PREVP
      SSFLAG         = SSFLAGP
C
C                         Set the scheme flag to zero
      SCHFLG = 0
C
C
 5    CONTINUE
C
C
C
C
C                         Set T equal to the input variable (the entropy
C                         and internal energy options are not implemented
C                         in this version)
      T = INPVAR(1)
      NSUBI = INPVAR(2)
      ETA_PO = INPVAR(3)
      ETA_NO = INPVAR(4)
C
C
C                         Calc the quantum concentration of nucleons
      NQ = 2.36D-4*T**1.5
C
C                         Calc the Fermi integral coefficent
      UQ = 20.721
C
      MQ = (T/UQ)**1.5
C
      KQ = ((T/UQ)**2.5)/(2.0*PI**2)
C
      LQ = UQ*(MQ**OVR53)/(3.0*(PI**2))
C
      ETAMAX = 0.95*FINV12(2.0*(PI**2)*BRYDNS/MQ)
C
      IF(ETA_PO.GE.ETAMAX) ETA_PO = ETAMAX-0.1
      IF(ETA_NO.GE.ETAMAX) ETA_NO = ETAMAX-0.1
      NI_MIN = DMAX1(4.5D-2,BRYDNS)
      IF(NSUBI.LT.NI_MIN) NSUBI = NI_MIN+1.0D-3
C
      TCFLAG = 0
C
      cflag = 0
C
      NEWFLG = 1
C
C                    Start Newton-Raphson iteration here
C
C
      DO 30 I=1,MAXIT,1
C
        IT_NUM = I
C                       Set the "Negative" flag
        NGFLAG = 0
C
C
C
        NNOUT = MQ*F_1_2(ETA_NO)/(2.0*PI**2)
        NPOUT = MQ*F_1_2(ETA_PO)/(2.0*PI**2)
C
        NOUT = NPOUT+NNOUT
C
c20        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)
        VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
C
c20        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+
c20     1    CC*(1.0+DD)*NOUT**DD+DELTAM)
        VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
        F32_NO = F_3_2(ETA_NO)
C
        F32_PO = F_3_2(ETA_PO)
C
c20        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(
c20     1    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*PV_PR(NPOUT,NNOUT)
C
        MUN_O = T*ETA_NO+VNOUT
C
        MUP_O = T*ETA_PO+VPOUT
C
        MUALFA = 2.0*MUN_O+2.0*MUP_O+BALPHA-BPROUT*V_ALFA
C
        IF(ABS(MUALFA/T).LT.30.0) THEN
          ALFDNS = 8.0*NQ*DEXP(MUALFA/T)
        ELSEIF((MUALFA/T).LT.-30.0) THEN
          ALFDNS = 0.0
        ELSE
          ALFDNS = 8.0*NQ*DEXP(3.0D1)
        ENDIF
C
C
C                   These statements take out the alfas if the
C                   alpha particle enable flag is not set
        IF(ALFLAG.NE.1) THEN
          ALFDNS = 0.0
          MUALFA = -300.0
        ENDIF
C
C
C
        EXALFA = 1.0-ALFDNS*V_ALFA
C
C
        BPRALF = ALFDNS*T
C
c---------------------------------------------------
C
C
C             Calculate fraction of space occupied by nuclei
        U_NUC = (BRYDNS-EXALFA*NOUT-4.0*ALFDNS)/
     1        (NSUBI-EXALFA*NOUT-4.0*ALFDNS)
C
C
C            Is volume occupied by nuclei within acceptable limits?
cc        IF((U_NUC.LT.0.0).OR.((U_NUC-1.0).GT.-1.0E-20)) THEN
cc        IF((U_NUC.LT.0.0).OR.(U_NUC.GT.0.996)) THEN
        IF((U_NUC.LT.1.0d-17).OR.(U_NUC.GT.0.996)) THEN
          NGFLAG = 1
          GOTO 29
        ENDIF
C
C
C            Volume exclusion factor due to nuclei
        EXCLU = 1.0-U_NUC
C
C
C            If calculated nucleon and alfa densities are larger
C            than the baryon density then reduce the eta's
        IF((EXCLU*EXALFA*NOUT+EXCLU*4.0*ALFDNS).GT.BRYDNS) THEN
          NGFLAG = 1
          GOTO 29
        ENDIF
C
C
C            Calculate the internal (inside nuclei) proton fraction
C
        X = (BRYDNS*YE-(1.0-U_NUC)*(EXALFA*NPOUT+2.0*ALFDNS))/
     1    (U_NUC*NSUBI)
        COMPX = 1.0-X
C
C
C            Is X within reasonable (but not necessarily correct)
C            limits? (YE may not be the lower bound on X !!!)
cccc        X_MIN = DMAX1(1.0D-2,(YE-0.05))
        X_MIN = DMAX1(1.0D-2,(0.8*YE))
cc        x_min = 0.95*ye
        IF((X.LT.X_MIN).OR.(X.GT.0.6)) THEN
          NGFLAG = 1
          GOTO 29
        ENDIF
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C                     Calculate critical temperature & its X derivative
        TSC_12 = 87.76*((COMP/375.0)**0.5)*((0.155/NSUBS)**OVR3)
c
cccdebug      tsc_12 = 1.0d8
c
        TSUBC = TSC_12*X*COMPX
        DTCDX = TSC_12*(1.0-2.0*X)
        DTCDXX = -2.0*TSC_12
        H = 1.0-2.0*(T/TSUBC)**2+(T/TSUBC)**4
C
cc        tsubc = tsc_12*0.25
cc        dtcdx = 0.0
cc        dtcdxx = 0.0
C
C
CC        TSUBC = 80.0*X*COMPX
CC        DTCDX = 80.0*(1.0-2.0*X)
CC        DTCDXX = -160.0
C
C                     If the X is such that T is greater than the
C                     critical temperature then fix NSUBI so that
C                     it lies in the bounds of acceptable parameter
C                     space
        ftflag = 0
        IF(((T.GT.TSUBC).OR.(H.LE.0.0)).AND.(SCHFLG.EQ.0)) THEN
C                       If this is an initial guess, then lower
C                       NSUBI untill we get a good X
          IF(NEWFLG.EQ.1) THEN
cc        write(*,*) ' nuc exceeded Tc'
cc        write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc
            ZNG = 2.0*(PI**2)*BRYDNS*(1.0-0.1*YE)/(1.5*MQ)
            ZPG = 2.0*(PI**2)*BRYDNS*0.1*YE/(1.5*MQ)
            IF(TCFLAG.NE.1) THEN
              TCFLAG = 1
              ETA_PO = FINV12(ZPG)-0.0
              ETA_NO = FINV12(ZNG)-0.0
            ELSE
              ETA_PO = ETA_PO-2.0/T
              ETA_NO = ETA_NO-2.0/T
              NSUBI = DMAX1(0.9*NSUBI,5.1D-2)
            ENDIF
            IF(DBFLAG.EQ.1) THEN
              WRITE(*,2000) '1',i,NSUBI,ETA_PO,ETA_NO,DNSUBI
            ENDIF
            GOTO 30
          ELSE
C                       Otherwise go back and cut the stepsize in
C                       half since it was obviously too big
            NGFLAG = 1
            GOTO 29
          ENDIF
        ELSEIF(((T.GT.TSUBC).OR.(H.LE.0.0)).AND.(SCHFLG.EQ.1)) THEN
          ftflag = 1
          tsubc = 80.0*(0.25+0.5*ye)*(0.75-0.25*ye)
C
        ENDIF
C
C
        R_0 = (0.75/(PI*NSUBS))**OVR3
        Q = (384.0*PI*(R_0**2)*SIG_S/SYM_S)-16.0
C
C                        Calculate surface functions of the internal
C                        (nuclear) proton fraction, X
        SIGMA = 1.0/(Q+1.0/(X**3)+1.0/(COMPX**3))
        OVRX4 = (1.0/X**4)-(1.0/COMPX**4)
        DSIGDX = 3.0*(SIGMA**2)*OVRX4
        SIGSGP = DSIGDX/SIGMA
        SIGSG2 = 18.0*(SIGMA**2)*OVRX4**2-12.0*SIGMA*((1.0/X**5)+
     1  (1.0/COMPX**5))
C
C                        If T is less than critical temp then
        IF(T.LT.TSUBC) THEN
C                        Calculate the surface energy temperature factor
C                        and its X and T derivatives
          H = 1.0-2.0*(T/TSUBC)**2+(T/TSUBC)**4
          HPRIM = -4.0*T/(TSUBC**2)+4.0*((T/TSUBC)**3)/TSUBC
          HPPRIM = -4.0/(TSUBC**2)+12.0*(T**2)/(TSUBC**4)
          DHDX = 4.0*(T**2/TSUBC**3-T**4/TSUBC**5)*DTCDX
          DHDXX = 4.0*(T**2/TSUBC**3-T**4/TSUBC**5)*DTCDXX+
     1    4.0*(-3.0*T**2/TSUBC**4+5.0*T**4/TSUBC**6)*(DTCDX**2)
          HX = DHDX/H
          DHDTDX = 8.0*(T/TSUBC**3-2.0*(T**3)/TSUBC**5)*DTCDX
C
C
C                        X independent version of TZERO
c          TZERO = 0.25*TSC_12
c          DTZDX = 0.0
c          DTZDXX = 0.0
C                        X dependent version of TZERO
c          TZERO = TSUBC
c          DTZDX = DTCDX
c          DTZDXX = DTCDXX
C
C
C
C                        Coulomb liquid correction factors and their
C                        derivatives
c          W = 1-(T/TZERO)**2
c          DWDX = 2.0*(T**2)*DTZDX/(TZERO**3)
c          DWDT = -2.0*T/(TZERO**2)
c          DWDTDX = 4.0*T*DTZDX/(TZERO**3)
c          DWDXDX = 2.0*(T**2)*
c     1    (DTZDXX/(TZERO**3)-3.0*(DTZDX**2)/(TZERO**4))
c          DWDTDT = -2.0/(TZERO**2)
C
          w = 1.0
          dwdt = 0.0
          dwdx = 0.0
          dwdtdx = 0.0
          dwdxdx = 0.0
          dwdtdt = 0.0
C
C
C
C                        Calc lattice factor & derivatives & products
C
          EXCLU = 1.0-U_NUC
          COMPU = 1.0-U_NUC
C
          DU = DMAX1(1.0D-15, (1.0-1.5*W*U_NUC**OVR3+0.5*U_NUC))
          DMU = DMAX1(1.0D-15,(1.0-1.5*W*(1.0-U_NUC+1.0E-20)**OVR3+
     1    0.5*(1.0-U_NUC)))
C
          DUP = -0.5*W*U_NUC**M2OVR3+0.5
          DMUP =-0.5*W*(1.0-U_NUC+1.0E-20)**M2OVR3+0.5
          DUPP = OVR3*W*((U_NUC+1.0D-20)**M5OVR3)
          DMUPP = OVR3*W*((1.0-U_NUC)+1.0E-20)**M5OVR3
C
C                Derivatives w.r.t. T
C
          DUT = -1.5*DWDT*U_NUC**OVR3
          DMUT = -1.5*DWDT*(1.0-U_NUC+1.0E-20)**OVR3
          DUPT = -0.5*DWDT*U_NUC**M2OVR3
          DMUPT = -0.5*DWDT*(1.0-U_NUC+1.0E-20)**M2OVR3
C
C                Derivatives w.r.t. X
C
          DUX = -1.5*DWDX*U_NUC**OVR3
          DMUX = -1.5*DWDX*(1.0-U_NUC+1.0E-20)**OVR3
          DUPX = -0.5*DWDX*U_NUC**M2OVR3
          DMUPX = -0.5*DWDX*(1.0-U_NUC+1.0E-20)**M2OVR3
C
C                Second derivatives w.r.t. X
C
          DUXX = -1.5*DWDXDX*U_NUC**OVR3
          DMUXX = -1.5*DWDXDX*(1.0-U_NUC+1.0E-20)**OVR3
C
C                Second derivatives w.r.t. T
C
          DUTT = -1.5*DWDTDT*U_NUC**OVR3
          DMUTT = -1.5*DWDTDT*(1.0-U_NUC+1.0E-20)**OVR3
C
C                Second derivatives w.r.t. X & T
C
          DUXT = -1.5*DWDTDX*U_NUC**OVR3
          DMUXT = -1.5*DWDTDX*(1.0-U_NUC+1.0E-20)**OVR3
C
C
          TMP1 = (U_NUC**2)+(COMPU**2)+0.6*(U_NUC*COMPU)**2
          TMP1P = 4.0*U_NUC-2.0+
     1    2.0*0.6*(U_NUC*COMPU**2-COMPU*U_NUC**2)
          TMP1PP = 4.0+2.0*0.6*(COMPU**2-4.0*U_NUC*COMPU+U_NUC**2)
C
          TMP2 = COMPU*(DU**OVR3)
          TMP2P = -1.0*DU**OVR3+OVR3*COMPU*(DU**M2OVR3)*DUP
          TMP2PP = -OVR23*(DU**M2OVR3)*DUP-OVR29*COMPU*
     1    (DU**M5OVR3)*DUP**2+OVR3*COMPU*(DU**M2OVR3)*DUPP
C
          TMP2T = OVR3*COMPU*(DU**M2OVR3)*DUT
          TMP2X = OVR3*COMPU*(DU**M2OVR3)*DUX
          TMP2XX = OVR3*COMPU*(DU**M2OVR3)*DUXX+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*(DUX**2)
          TMP2TT = OVR3*COMPU*(DU**M2OVR3)*DUTT+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*(DUT**2)
          TMP2XT = OVR3*COMPU*(DU**M2OVR3)*DUXT+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUX*DUT
          TMP2PT = -OVR3*(DU**M2OVR3)*DUT+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUP*DUT+
     2        OVR3*COMPU*(DU**M2OVR3)*DUPT
          TMP2PX = -OVR3*(DU**M2OVR3)*DUX+
     1        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUP*DUX+
     2        OVR3*COMPU*(DU**M2OVR3)*DUPX
C
C
C
          TMP3 = U_NUC*(DMU**OVR3)
          TMP3P = (DMU**OVR3)-OVR3*U_NUC*(DMU**M2OVR3)*DMUP
          TMP3PP = -OVR23*(DMU**M2OVR3)*DMUP-OVR29*U_NUC*
     1    (DMU**M5OVR3)*(DMUP**2)+OVR3*U_NUC*(DMU**M2OVR3)*DMUPP
C
          TMP3T = OVR3*U_NUC*(DMU**M2OVR3)*DMUT
          TMP3X = OVR3*U_NUC*(DMU**M2OVR3)*DMUX
          TMP3XX = OVR3*U_NUC*(DMU**M2OVR3)*DMUXX+
     1        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*(DMUX**2)
          TMP3TT = OVR3*U_NUC*(DMU**M2OVR3)*DMUTT+
     1        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*(DMUT**2)
          TMP3XT = OVR3*U_NUC*(DMU**M2OVR3)*DMUXT+
     1        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*DMUX*DMUT
          TMP3PT = OVR3*(DMU**M2OVR3)*DMUT-OVR3*M2OVR3*U_NUC*
     1      (DMU**M5OVR3)*DMUP*DMUT-OVR3*U_NUC*(DMU**M2OVR3)*DMUPT
 
          TMP3PX = OVR3*(DMU**M2OVR3)*DMUX-OVR3*M2OVR3*U_NUC*
     1      (DMU**M5OVR3)*DMUP*DMUX-OVR3*U_NUC*(DMU**M2OVR3)*DMUPX
C
C
C                 Combination D function
C
          SCRDU = U_NUC*COMPU*(TMP2+TMP3)/TMP1
          SCRDUT = U_NUC*COMPU*(TMP2T+TMP3T)/TMP1
          SCRDUX = U_NUC*COMPU*(TMP2X+TMP3X)/TMP1
          SCRDXX = U_NUC*COMPU*(TMP2XX+TMP3XX)/TMP1
          SCRDTT = U_NUC*COMPU*(TMP2TT+TMP3TT)/TMP1
          SCRDXT = U_NUC*COMPU*(TMP2XT+TMP3XT)/TMP1
C
          SCRD = SCRDU/U_NUC
          SCRDT = SCRDUT/U_NUC
          SCRDX = SCRDUX/U_NUC
C
          SCRD2 = SCRDU/COMPU
          SCRD2T = SCRDUT/COMPU
          SCRD2X = SCRDUX/COMPU
C
          SCRDUP = SCRD-SCRD2+U_NUC*COMPU*
     1    ((TMP2P+TMP3P)/TMP1-(TMP2+TMP3)*TMP1P/TMP1**2)
C
          SCRDPT = SCRDT-SCRD2T+U_NUC*COMPU*
     1    ((TMP2PT+TMP3PT)/TMP1-(TMP2T+TMP3T)*TMP1P/TMP1**2)
C
          SCRDPX = SCRDX-SCRD2X+U_NUC*COMPU*
     1    ((TMP2PX+TMP3PX)/TMP1-(TMP2X+TMP3X)*TMP1P/TMP1**2)
C
          SCRDPP = (SCRDUP-SCRD)/U_NUC-(SCRD2+SCRDUP)/COMPU+
     1    (1.0-2.0*U_NUC)*
     2    ((TMP2P+TMP3P)/TMP1-(TMP2+TMP3)*TMP1P/TMP1**2)+U_NUC*COMPU*
     3    ((TMP2PP+TMP3PP)/TMP1-2.0*(TMP2P+TMP3P)*TMP1P/TMP1**2-
     4    (TMP2+TMP3)*TMP1PP/TMP1**2+
     5    2.0*(TMP2+TMP3)*(TMP1P**2)/TMP1**3)
C
C
c
c           bubble D function
cbub          scrdu = (1.0-u_nuc)*dmu**ovr3
cbub          scrd = scrdu/u_nuc
cbub          scrd2 = dmu**ovr3
cbub          scrdup = -1.0*dmu**ovr3-
cbub     1    ovr3*(1.0-u_nuc)*dmup*dmu**m2ovr3
cbub          scrdpp = ovr23*dmup*dmu**m2ovr3-ovr29*(1.0-u_nuc)*
cbub     1    dmu**m5ovr3*dmup**2+ovr3*(1.0-u_nuc)*dmu**m2ovr3*dmupp
c
c
c           nuclei D function
cnuc          scrdu = u_nuc*du**ovr3
cnuc          scrd = du**ovr3
cnuc          scrd2 = scrdu/(1.0-u_nuc)
cnuc          scrdup = du**ovr3+ovr3*u_nuc*dup*du**m2ovr3
cnuc          scrdpp = ovr23*dup*du**m2ovr3-ovr29*u_nuc*
cnuc     1    (du**m5ovr3)*(dup**2)+ovr3*u_nuc*(du**m2ovr3)*dupp
c
c
C
          ZETA_0 = CSSCAL*6.035204*(SIG_S*(16.0+Q))**OVR23
C
C                        Surface energy coefficent
          ZETA = ZETA_0*(H*SIGMA*X*NSUBI)**OVR23
C
C                        Derivative of Zeta w.r.t. X
          DZDT = OVR23*ZETA*HPRIM/H
C
C                        Derivative of Zeta w.r.t. X
          DZDX = OVR23*ZETA*(DHDX/H+SIGSGP+1.0/X)
C
C                        Derivative of Zeta w.r.t. NSUBI
          DZDNI = OVR23*ZETA/NSUBI
C
C
C
C                        Nuclear radius
          RSUBN = 9.0*H*SIGMA*SIG_S*(16.0D0+Q)*U_NUC*(1.0-U_NUC)/
     1    (2.0*ZETA*SCRDU)
C
C                        Nuclear volume
          VSUBN = 4.0*PI*(RSUBN**3)/3.0
C
C                        Atomic number
          A = NSUBI*VSUBN
C
C                        Now calc surface, Coulomb free energies
C
          FSUBSC = ZETA*SCRDU/BRYDNS
          FSUBS = OVR23*ZETA*SCRDU/BRYDNS
          FSUBC = OVR3*ZETA*SCRDU/BRYDNS
C
C
C
C                   Translational chemical potential
          MUSUBT = TRSCAL*
     1        T*DLOG((1.0-U_NUC)*(U_NUC*NSUBI)/(NQ*AZERO**2.5))
C
C                   Derivative of trans. chem. potential w.r.t. T
          DMUTDT = TRSCAL*(MUSUBT/T-1.5)
C
C                   Translational free energy per baryon
          FTRANS = TRSCAL*H*(MUSUBT-T)/AZERO
C
C                            if T is above the critical temperature
        ELSE
          A = 0.0
          RSUBN = 0.0
          VSUBN = 0.0
          FSUBS = 0.0
          FSUBC = 0.0
          FTRANS = 0.0
        ENDIF
C                            Calc ratio of NSUBI to NSUBS
        NRATIO = NSUBI/NSUBS
C
C
c20        VNI = 2.0*AA*NSUBI+4.0*BB*X*NSUBI+CC*(1.0+DD)*NSUBI**DD
        VNI = PVN(X*NSUBI,(1.0-X)*NSUBI)
C
c20        VPI = 2.0*AA*NSUBI+4.0*BB*(1.0-X)*NSUBI+
c20     1    CC*(1.0+DD)*NSUBI**DD+DELTAM
        VPI = PVP(X*NSUBI,(1.0-X)*NSUBI)
C
c---------------------------------------------------
C
        ZNI = 2.0*(PI**2)*NSUBI*(1.0-X)/MQ
C
        ZPI = 2.0*(PI**2)*NSUBI*X/MQ
C
        ETA_NI = FINV12(ZNI)
C
        ETA_PI = FINV12(ZPI)
C
        MUN_I = T*ETA_NI+VNI
C
        MUP_I = T*ETA_PI+VPI
C
        F32_NI = F_3_2(ETA_NI)
C
        F32_PI = F_3_2(ETA_PI)
C
c20        PSUBI = LQ*(F32_NI+F32_PI)+
c20     1    (NSUBI**2)*(AA+4.0*BB*X*(1.0-X))+DD*CC*NSUBI**(1.0+DD)
        PSUBI = LQ*(F32_NI+F32_PI)+PV_PR(X*NSUBI,(1.0-X)*NSUBI)
C
C
        BN = OVR23*ZETA*SCRD*(SIGSGP+HX+1.5*SCRDUX/SCRDU)*X/NSUBI-
     1  TRSCAL*(1.0-U_NUC)*(MUSUBT*(H-X*DHDX)/AZERO+X*DHDX*T/AZERO)
C
        BP = -OVR23*ZETA*SCRD*
     1 ((SIGSGP+HX+1.5*SCRDUX/SCRDU)*COMPX+1.0/X)/NSUBI-
     1 TRSCAL*(1.0-U_NUC)*
     2 (MUSUBT*(H+DHDX*COMPX)/AZERO-DHDX*T*COMPX/AZERO)
C
        BSUBP = ZETA*SCRDUP-OVR23*ZETA*SCRD-
     1        TRSCAL*U_NUC*NSUBI*H*MUSUBT/AZERO
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
cc        GPI = 2.0*FHALFI(ETA_PI)
cc        GPI = 2.0*FHALFI(2.0*(pi**2)*x*nsubi/mq)
cc        GNI = 2.0*FHALFI(ETA_NI)
cc        GNI = 2.0*FHALFI(2.0*(pi**2)*(1.0-x)*nsubi/mq)
c
cc        GPO = 2.0*FHALFO(ETA_PO)
cc        GNO = 2.0*FHALFO(ETA_NO)
C
c
        GPO = 2.0*FHALF(ETA_PO)
        GNO = 2.0*FHALF(ETA_NO)
        GPI = 2.0*FHALF(ETA_PI)
        GNI = 2.0*FHALF(ETA_NI)
C
C                  Derivatives of inside potentials
C
c20        DVPIDP = 2.0*AA+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0))
        DVPIDP = DPVPDP(X*NSUBI,(1.0-X)*NSUBI)
c20        DVPIDN = 2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0))
        DVPIDN = DPVPDN(X*NSUBI,(1.0-X)*NSUBI)
c20        DVNIDP = DVPIDN
        DVNIDP = DPVNDP(X*NSUBI,(1.0-X)*NSUBI)
c20        DVNIDN = DVPIDP
        DVNIDN = DPVNDN(X*NSUBI,(1.0-X)*NSUBI)
C
C                  Derivatives of outside potentials
C
c20        DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) )
        DVPODP = EIFLAG*DPVPDP(NPOUT,NNOUT)
c20        DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)))
        DVPODN = EIFLAG*DPVPDN(NPOUT,NNOUT)
c20        DVNODP = DVPODN
        DVNODP = EIFLAG*DPVNDP(NPOUT,NNOUT)
c20        DVNODN = DVPODP
        DVNODN = EIFLAG*DPVNDN(NPOUT,NNOUT)
C
C                  Derivatives of inside K.E. densities
C
        MSSCON = 3.0*MASSN/((HBAR*C)**2)
        DTPIDP = MSSCON*T*GPI
        DTPIDN = 0.0
        DTNIDP = 0.0
        DTNIDN = MSSCON*T*GNI
C
C                  Derivatives of outside K.E. densities
C
        DTPODP = MSSCON*T*GPO
        DTPODN = 0.0
        DTNODP = 0.0
        DTNODN = MSSCON*T*GNO
C
C
C                  Derivatives of inside chem. potentials
C
        DMPIDP = T*GPI/(X*NSUBI)+DVPIDP
        DMPIDN = DVPIDN
        DMNIDP = DVNIDP
        DMNIDN = T*GNI/((1.0-X)*NSUBI)+DVNIDN
C
C                  Derivatives of outside chem. potentials
C
        DMPODP = T+DVPODP*NPOUT/GPO
        DMPODN = DVPODN*NNOUT/GNO
        DMNODP = DVNODP*NPOUT/GPO
        DMNODN = T+DVNODN*NNOUT/GNO
C
C                  Derivatives of inside pressure
C
        DPIDP = X*NSUBI*DMPIDP+(1.0-X)*NSUBI*DMNIDP
        DPIDN = X*NSUBI*DMPIDN+(1.0-X)*NSUBI*DMNIDN
C
C                  Derivatives of outside pressure
C
        DPODP = NPOUT*DMPODP+NNOUT*DMNODP
        DPODN = NPOUT*DMPODN+NNOUT*DMNODN
C
C                  Derivatives of alpha pressure
C
        DPADP = ALFDNS*
     1  ( (2.0-NPOUT*V_ALFA)*DMPODP+(2.0-NNOUT*V_ALFA)*DMNODP )
        DPADN = ALFDNS*
     1  ( (2.0-NPOUT*V_ALFA)*DMPODN+(2.0-NNOUT*V_ALFA)*DMNODN )
C
C
        N1 = NSUBI-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS
        N2 = NSUBI*X-EXALFA*NPOUT-2.0*ALFDNS
C
C                  Derivatives of U
C
        DUDPO = -EXCLU*(EXALFA*NPOUT/GPO+
     1           (4.0-NOUT*V_ALFA)*DPADP/T)/N1
        DUDNO = -EXCLU*(EXALFA*NNOUT/GNO+
     1           (4.0-NOUT*V_ALFA)*DPADN/T)/N1
        DUDNI = -U_NUC/N1
C
C                  Derivatives of X
C
        DXDPO = -(N2*DUDPO+EXCLU*(EXALFA*NPOUT/GPO+
     1           (2.0-NPOUT*V_ALFA)*DPADP/T))/(U_NUC*NSUBI)
        DXDNO = -(N2*DUDNO+EXCLU*(2.0-NPOUT*V_ALFA)*DPADN/T)/
     1           (U_NUC*NSUBI)
        DXDNI = (N2-X*N1)/(NSUBI*N1)
C
C                  Derivatives of B's w.r.t. NSUBI
C
        DB1DNI = TRSCAL*( -U_NUC*H*(MUSUBT+T)/AZERO )+
     1      OVR23*ZETA*(SCRDUP-OVR23*SCRD)/NSUBI
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X
C
        DB2DNI = -2.0*ZETA*SCRD*TMP4/(9.0*NSUBI**2)-
     1  TRSCAL*( (COMPU*T/(AZERO*NSUBI))*(H+COMPX*DHDX) )
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)
        DB3DNI = -2.0*ZETA*SCRD*X*TMP4/(9.0*NSUBI**2)-
     1          TRSCAL*( ((COMPU*T)/(AZERO*NSUBI))*(H-X*DHDX) )
C
c
c
C                  Derivatives of B's w.r.t. X
C
        DB1DX = OVR23*ZETA*(SCRDUP-OVR23*SCRD)*(SIGSGP+DHDX/H+1.0/X)+
     1  OVR23*ZETA*(SCRDPX-OVR23*SCRDX)-
     2  TRSCAL*( U_NUC*NSUBI*DHDX*MUSUBT/AZERO )
C
C
C
C
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X
C
        TMP5 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU+(X**(-2))+(X-1.0)*
     1  (SIGSG2-(SIGSGP**2)-(DHDX/H)**2+DHDXX/H+1.5*SCRDXX/SCRDU-
     2  1.5*(SCRDUX/SCRDU)**2)
C
C
        DB2DX = OVR23*(ZETA*SCRDUX+SCRDU*DZDX)*TMP4/(U_NUC*NSUBI)+
     1      OVR23*ZETA*SCRD*TMP5/NSUBI-TRSCAL*(
     2      COMPU*(DHDX*MUSUBT+(DHDXX*(1.0-X)-DHDX)*(MUSUBT-T))/AZERO)
C
C
C
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*X
C
        TMP5 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU+X*
     1         (SIGSG2-(SIGSGP**2)-(DHDX/H)**2+DHDXX/H+
     2       1.5*SCRDXX/SCRDU-1.5*(SCRDUX/SCRDU)**2)
C
        DB3DX = OVR23*(ZETA*SCRDUX+SCRDU*DZDX)*TMP4/(U_NUC*NSUBI)+
     1      OVR23*ZETA*SCRD*TMP5/NSUBI-
     2      TRSCAL*( COMPU*(DHDX*T-X*DHDXX*(MUSUBT-T))/AZERO )
C
C
C
C                  Derivatives of B's w.r.t. U_NUC
C
        DB1DU = ZETA*(SCRDPP-OVR23*SCRDUP/U_NUC+OVR23*SCRD/U_NUC)-
     1  TRSCAL*( NSUBI*H*(MUSUBT+T*(1.0-2.0*U_NUC)/(1.0-U_NUC))/AZERO )
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X
        TMP5 = (X-1.0)*1.5*(SCRDPX/SCRDU-SCRDUX*SCRDUP/SCRDU**2)
        DB2DU = (OVR23*ZETA*SCRD/NSUBI)*TMP4*(SCRDUP/SCRDU-1.0/U_NUC)+
     1    OVR23*ZETA*SCRDU*TMP5/(U_NUC*NSUBI)+
     1    TRSCAL*( (H*MUSUBT+DHDX*COMPX*(MUSUBT-T))/AZERO-
     2    (T*(1.0-2.0*U_NUC)/U_NUC)*(H+DHDX*COMPX)/AZERO )
C
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*X
        TMP5 = X*1.5*(SCRDPX/SCRDU-SCRDUP*SCRDUX/SCRDU**2)
        DB3DU = OVR23*ZETA*SCRD*TMP4*(U_NUC*SCRDUP/SCRDU-1.0)/
     1 (U_NUC*NSUBI)+OVR23*ZETA*SCRDU*TMP5/(U_NUC*NSUBI)+
     2  TRSCAL*( (H*MUSUBT-X*DHDX*(MUSUBT-T))/AZERO-
     3 T*(1.0-2.0*U_NUC)*(H-X*DHDX)/(AZERO*U_NUC) )
C
C
C                      A1 derivatives
C
        DA1ID1 = X*DPIDP+(1.0-X)*DPIDN+NSUBI*(DPIDP-DPIDN)*DXDNI
        DA1ID2 = NSUBI*(DPIDP-DPIDN)*DXDPO
        DA1ID3 = NSUBI*(DPIDP-DPIDN)*DXDNO
C
        DA1OD1 = 0.0
        DA1OD2 = DPODP+DPADP
        DA1OD3 = DPODN+DPADN
C
        DB1D1 = DB1DNI+DB1DX*DXDNI+DB1DU*DUDNI
        DB1D2 = DB1DX*DXDPO+DB1DU*DUDPO
        DB1D3 = DB1DX*DXDNO+DB1DU*DUDNO
C
        DA1D1 = DA1ID1-DB1D1-DA1OD1
        DA1D2 = DA1ID2-DB1D2-DA1OD2
        DA1D3 = DA1ID3-DB1D3-DA1OD3
C
C                      A3 derivatives
C
        DA3ID1 = X*DMNIDP+(1.0-X)*DMNIDN+NSUBI*(DMNIDP-DMNIDN)*DXDNI
        DA3ID2 = NSUBI*(DMNIDP-DMNIDN)*DXDPO
        DA3ID3 = NSUBI*(DMNIDP-DMNIDN)*DXDNO
C
        DA3OD1 = 0.0
        DA3OD2 = DMNODP
        DA3OD3 = DMNODN
C
        DB3D1 = DB3DNI+DB3DX*DXDNI+DB3DU*DUDNI
        DB3D2 = DB3DX*DXDPO+DB3DU*DUDPO
        DB3D3 = DB3DX*DXDNO+DB3DU*DUDNO
C
        DA3D1 = DA3ID1-DB3D1-DA3OD1
        DA3D2 = DA3ID2-DB3D2-DA3OD2
        DA3D3 = DA3ID3-DB3D3-DA3OD3
C
C                      A2 derivatives
C
        DA2ID1 = X*DMPIDP+(1.0-X)*DMPIDN+NSUBI*(DMPIDP-DMPIDN)*DXDNI
        DA2ID2 = NSUBI*(DMPIDP-DMPIDN)*DXDPO
        DA2ID3 = NSUBI*(DMPIDP-DMPIDN)*DXDNO
C
        DA2OD1 = 0.0
        DA2OD2 = DMPODP
        DA2OD3 = DMPODN
C
        DB2D1 = DB2DNI+DB2DX*DXDNI+DB2DU*DUDNI
        DB2D2 = DB2DX*DXDPO+DB2DU*DUDPO
        DB2D3 = DB2DX*DXDNO+DB2DU*DUDNO
C
        DA2D1 = DA2ID1-DB2D1-DA2OD1
        DA2D2 = DA2ID2-DB2D2-DA2OD2
        DA2D3 = DA2ID3-DB2D3-DA2OD3
C
C
C                      Eta derivatives
C
        DNDETN = NNOUT/GNO
        DPDETP = NPOUT/GPO
C
        DA1DN = DA1D1
        DA1ETP = DA1D2
        DA1ETN = DA1D3
C
        DA2DN = DA2D1
        DA2ETP = DA2D2
        DA2ETN = DA2D3
C
        DA3DN = DA3D1
        DA3ETP = DA3D2
        DA3ETN = DA3D3
C
C
C
C
        A1 = PSUBI-BSUBP-BPROUT-BPRALF
        A2 = MUP_I-BP-MUP_O
        A3 = MUN_I-BN-MUN_O
C
C
C                          Unset the "new" flag
        NEWFLG = 0
C
        DETERM = DA1DN*(DA2ETP*DA3ETN-DA2ETN*DA3ETP)-
     1           DA1ETP*(DA2DN*DA3ETN-DA2ETN*DA3DN)+
     2           DA1ETN*(DA2DN*DA3ETP-DA2ETP*DA3DN)
C
        DNSUBI = -1.0*(A1*(DA2ETP*DA3ETN-DA2ETN*DA3ETP)+
     1           A2*(DA3ETP*DA1ETN-DA1ETP*DA3ETN)+
     2           A3*(DA1ETP*DA2ETN-DA1ETN*DA2ETP))/DETERM
C
C
        DETAP = -1.0*(A1*(DA2ETN*DA3DN-DA2DN*DA3ETN)+
     1          A2*(DA1DN*DA3ETN-DA1ETN*DA3DN)+
     2          A3*(DA1ETN*DA2DN-DA1DN*DA2ETN))/DETERM
C
C
        DETAN = -1.0*(A1*(DA2DN*DA3ETP-DA2ETP*DA3DN)+
     1          A2*(DA1ETP*DA3DN-DA1DN*DA3ETP)+
     2          A3*(DA1DN*DA2ETP-DA1ETP*DA2DN))/DETERM
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C
C                        Check the step size in NSUBI
        IF(ABS(DNSUBI/NSUBI).GT.0.04) THEN
          DNSUBI = 0.04*DNSUBI*NSUBI/ABS(DNSUBI)
        ENDIF
 26     CONTINUE
        NSUBIN = NSUBI+DNSUBI
        IF((NSUBIN.LT.DMAX1(4.5D-2,BRYDNS)).OR.(NSUBIN.GT.0.25)) THEN
          DNSUBI = 0.5*DNSUBI
          GOTO 26
        ENDIF
C
C                        Check the step size in ETA_PO
        IF(ABS(DETAP).GT.4.0) THEN
          DETAP = 4.0*DETAP/ABS(DETAP)
        ENDIF
 27     CONTINUE
        NETAP = ETA_PO+DETAP
        IF((NETAP.LT.-5000.0).OR.(NETAP.GT.ETAMAX)) THEN
          DETAP = 0.5*DETAP
          GOTO 27
        ENDIF
C
C                        Check the step size in ETA_NO
        IF(ABS(DETAN).GT.4.0) THEN
          DETAN = 4.0*DETAN/ABS(DETAN)
        ENDIF
 28     CONTINUE
        NETAN = ETA_NO+DETAN
        IF((NETAN.LT.-5000.0).OR.(NETAN.GT.ETAMAX)) THEN
          DETAN = 0.5*DETAN
          GOTO 28
        ENDIF
C
C
C                        Update the variables
ccc        if(i.lt.30) write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc
 1205   format(i3,9(1pe21.14))
c
        NSUBI = NSUBI+DNSUBI
        ETA_PO = ETA_PO+DETAP
        ETA_NO = ETA_NO+DETAN
C
C
C
C                        If the required tolarences have been met
C                        break out of the loop
        IF((ABS(DNSUBI).LT.NSIACC).AND.(ABS(DETAP).LT.PRTACC)
     1    .AND.(ABS(DETAN).LT.NUTACC) ) THEN
          GOTO 40
        ELSE
      IF(DBFLAG.EQ.1) THEN
        WRITE(*,2000) '2',i,NSUBI,ETA_PO,ETA_NO,DNSUBI
      ENDIF
          GOTO 30
        ENDIF
C
C
 29     CONTINUE
        IF(NEWFLG.NE.1) THEN
          cflag = cflag+1
          DNSUBI = 0.5*DNSUBI
          NSUBI = NSUBI-DNSUBI
          DETAP = 0.5*DETAP
          ETA_PO = ETA_PO-DETAP
          DETAN = 0.5*DETAN
          ETA_NO = ETA_NO-DETAN
          IF(DBFLAG.EQ.1) THEN
            WRITE(*,2000) '3',i,NSUBI,ETA_PO,ETA_NO,DNSUBI
          ENDIF
          GOTO 30
        ELSE
          NSUBI = NSUBS
cc          ETA_PO = ETA_PO-0.5/T
cc          ETA_NO = ETA_NO-0.5/T
          ETA_PO = ETA_PO-2.0/T
          ETA_NO = ETA_NO-2.0/T
        ENDIF
C
C
      IF(DBFLAG.EQ.1) THEN
        WRITE(*,2000) '4',i,NSUBI,ETA_PO,ETA_NO,DNSUBI
      ENDIF
 2000   FORMAT(t2,a,1x,i3,1x,f8.5,3(1X,G13.5))
C
C
 30   CONTINUE
C
C            If scheme 1 has failed try scheme 2
      if(schflg.eq.0) then
        schflg = 1
        goto 5
      endif
c
c
      SSFLAG = 0
      GOTO 999
C
C                    Branch label to break out of DO 30 iteration
 40   CONTINUE
C
C
C                    The following logic determines whether this was
C                    the correct scheme to use, and if not then which
C                    one should be used
C
      if(ftflag.ne.0) then
        ssflag = 4
        goto 999
      endif
C
C                    If calculated critical temperature is less than T,
C                    then switch to the scheme with no nuclei
      IF(T.GE.TSUBC) THEN
C                    Set flag to indicate FAILURE
        SSFLAG = 0
        GOTO 999
      ENDIF
C
C
C                    If fraction of nuclei present is zero and no switch
C                    has been made then switch to the no nuclei scheme
      IF(U_NUC.LE.0.0) THEN
C                    Set flag to indicate FAILURE
        SSFLAG = 0
        GOTO 999
      ELSEIF(U_NUC.GT.1.0) THEN
C                    Set flag to indicate FAILURE
        SSFLAG = 0
        GOTO 999
      ELSE
C                    Set flag to indicate success
        SSFLAG = 1
      ENDIF
C
C
C
C                    If eqns aren't really zeroed then fail
C
C
      IF( (ABS(A1).GT.1.0D-5).OR.(ABS(A2).GT.1.0D-5).OR.
     1    (ABS(A3).GT.1.0D-5) ) THEN
        SSFLAG = 0
cc        WRITE(*,*) ' NUCEOS: False convg; A = ',A1,A2,A3
        GOTO 999
      ENDIF
C
C
C
C
      IF(NSUBI.LT.0.05) THEN
        WRITE(*,*) 'NUCEOS:: <<WARNING>> NSUBI GETTING CLOSE TO LB'
      ENDIF
C
C
C
      ZNI = 2.0*(PI**2)*NSUBI*(1.0-X)/MQ
C
      ZPI = 2.0*(PI**2)*NSUBI*X/MQ
C
      ETA_NI = FINV12(ZNI)
C
      ETA_PI = FINV12(ZPI)
C
      MUN_I = T*ETA_NI+VNI
C
      MUP_I = T*ETA_PI+VPI
C
      F32_NI = F_3_2(ETA_NI)
C
      F32_PI = F_3_2(ETA_PI)
C
      EXCLU = 1.0-U_NUC
      EXALFA = 1.0-ALFDNS*V_ALFA
C
C
C
C                    Calculate particle fractions
C
      XALFA = 4.0*EXCLU*ALFDNS/BRYDNS
      XNUT = NNOUT*EXCLU*EXALFA/BRYDNS
      XPROT = NPOUT*EXCLU*EXALFA/BRYDNS
      XH = 1.0-XPROT-XNUT-XALFA
      XHCHK = U_NUC*NSUBI/BRYDNS
C
      IF((XH.LT.HEAVCT).OR.(XHCHK.LT.HEAVCT)) THEN
C                    Set flag to indicate switch is being made
        SSFLAG = 0
cc        write(*,*) ' xh,xhchk = ',xh,xhchk
        GOTO 999
      ENDIF
C
      IF((XALFA.LT.0.0).OR.(XH.LT.0.0).OR.
     1   (XNUT.LT.0.0).OR.(XPROT.LT.0.0)) THEN
        SSFLAG = 0
        write(*,*) ' Xs hnpa = ',xh,xnut,xprot,xalfa
        GOTO 999
      ENDIF
C
C
C
C
C                    Baryons
C
C
      MUPROT = MUP_O
      MUN = MUN_O
      MUHAT = MUN-MUPROT
C
C
      IF(ABS((XH-XHCHK)/XHCHK).GT.1.0D-4) THEN
        SSFLAG = 0
        GOTO 999
CCC        WRITE(*,*) ' INCONSISTENCEY IN XH AT',T,BRYDNS,YE,XH,XHCHK
      ENDIF
C
      NUCDNS = BRYDNS*XH
C
      TAU_PO = KQ*F32_PO
      TAU_PI = KQ*F32_PI
C
      TAU_NO = KQ*F32_NO
      TAU_NI = KQ*F32_NI
C
      IF(NOUT.GT.0.0) XOUT = NPOUT/NOUT
C
C
C                    Calculate internal energy of outside nucleons,
C                    alpha particles, and nuclei (per baryon)
C
c20      BUOUT = (EXCLU*EXALFA/BRYDNS)*( UQ*(TAU_PO+TAU_NO)+EIFLAG*
c20     1    ( (NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+
c20     2    CC*NOUT**(1.0+DD)+NPOUT*DELTAM) )
      BUOUT = (EXCLU*EXALFA/BRYDNS)*(
     1    UQ*(TAU_PO+TAU_NO)+EIFLAG*PV_E(NPOUT,NNOUT) )
C
c20      BUNUC = XH*( ( UQ*(TAU_PI+TAU_NI)+(NSUBI**2)*
c20     1 (AA+4.0*BB*X*(1.0-X))+CC*NSUBI**(1.0+DD)+X*NSUBI*DELTAM )/
c20     2 NSUBI)+FSUBSC*(1.0-T*(SCRDUT/SCRDU+OVR23*HPRIM/H))+
c20     3 TRSCAL*
c20     4 (1.0-U_NUC)*XH*(FTRANS*(1.0-T*HPRIM/H)-H*(MUSUBT-2.5*T)/AZERO)
      BUNUC = XH*( (UQ*(TAU_PI+TAU_NI)+PV_E(X*NSUBI,(1.0-X)*NSUBI))/
     2 NSUBI)+FSUBSC*(1.0-T*(SCRDUT/SCRDU+OVR23*HPRIM/H))+
     3 TRSCAL*
     4 (1.0-U_NUC)*XH*(FTRANS*(1.0-T*HPRIM/H)-H*(MUSUBT-2.5*T)/AZERO)
C
C
      BUALFA = 0.25*XALFA*(1.5*T-BALPHA)
C
      BU = BUOUT+BUALFA+BUNUC
C
C
      BSOUT = (EXCLU*EXALFA/BRYDNS)*( (5.0*UQ/(3.0*T))*(TAU_NO+TAU_PO)-
     1 NNOUT*ETA_NO-NPOUT*ETA_PO )
C
C
C                    Calculate entropy of alpha particles (per baryon)
      BSALFA = -0.25*XALFA*(MUALFA/T-2.5)
C
C
      BSNUC = XH*( (5.0*UQ/(3.0*T))*(TAU_NI+TAU_PI)-
     1 NSUBI*(1.0-X)*ETA_NI-NSUBI*X*ETA_PI )/NSUBI-
     2 FSUBSC*(SCRDUT/SCRDU+OVR23*HPRIM/H)-
     3 XH*TRSCAL*(1.0-U_NUC)*
     4 ((FTRANS*HPRIM/H)+H*(MUSUBT/T-2.5)/AZERO)
C
C                    Calculate total baryon entropy (per baryon)
      BS = BSOUT+BSNUC+BSALFA
C
C                    Calculate free energy of outside nucleons (per baryon)
      BFOUT = BUOUT-T*BSOUT
C
C                    Calculate free energy of alpha particles (per baryon)
      BFALFA = BUALFA-T*BSALFA
C
C                    Calculate free energy of nuclei (per baryon)
      BFNUC = BUNUC-T*BSNUC
C
C                    Calculate total baryon free energy (per baryon)
      BFTOT = BFOUT+BFNUC+BFALFA
C
C                    Calculate pressure due to nuclei
      BPRNUC = -ZETA*(SCRDU-U_NUC*SCRDUP)+
     1 TRSCAL*U_NUC*NSUBI*H*((1.0-U_NUC)*T-U_NUC*MUSUBT)/AZERO
C
C
C                    Calculate total baryon pressure
      BPRESS = BPROUT+BPRALF+BPRNUC
C
C
C                    Leptons & Photons
C
      CALL EL_EOS(T,YE,BRYDNS)
C
C
C
C                    Total pressure and eng/ent per baryon
C
      FBARY = BFTOT+FSUBE
      PBARY = BPRESS+EPRESS
      MUBARY = YE*MUPROT+(1.0-YE)*MUN
      MU_MAT = YE*(MUPROT+MUSUBE)+(1.0-YE)*MUN
C
      FTOT = BFTOT+FSUBE+PF
      UTOT = BU+EU+PU
      STOT = BS+ES+PS
      PTOT = BPRESS+EPRESS+PPRESS
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------
C                Derivatives of thermodynamic variables
C-----------------------------------------------------------------------
C
C                 ------------------------------------
C                 !      Derivatives of exterior     !
C                 !      quantities                  !
C                 !      (w.r.t. Temp. and ETA's)    !
C                 !                                  !
C                 ------------------------------------
C
C
C                  Derivatives of exterior potentials
C                  w.r.t. particle densities
c20      DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) )
      DVPODP = EIFLAG*DPVPDP(NPOUT,NNOUT)
c20      DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)))
      DVPODN = EIFLAG*DPVPDN(NPOUT,NNOUT)
c20      DVNODP = DVPODN
      DVNODP = EIFLAG*DPVNDP(NPOUT,NNOUT)
c20      DVNODN = DVPODP
      DVNODN = EIFLAG*DPVNDN(NPOUT,NNOUT)
C
C
C                  Derviatives of exterior chem. pot. w.r.t. ETA's
C                  (at fixed T)
      DMPDEP = T+DVPODP*NPOUT/GPO
      DMPDEN = DVPODN*NNOUT/GNO
      DMNDEP = DVNODP*NPOUT/GPO
      DMNDEN = T+DVNODN*NNOUT/GNO
C
C                  Derivatives of pressure potential w.r.t.
C                  particle densities
c20      DV_DPO = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*DD*(1.0+DD)*(NOUT**DD) )
      DV_DPO = EIFLAG*DPVRDP(NPOUT,NNOUT)
c20      DV_DNO = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NPOUT+CC*DD*(1.0+DD)*(NOUT**DD) )
      DV_DNO = EIFLAG*DPVRDN(NPOUT,NNOUT)
C
C                  Derivatives of pressure potential w.r.t. ETA's
C                  (at fixed T)
      DV_DEP = DV_DPO*NPOUT/GPO
      DV_DEN = DV_DNO*NNOUT/GNO
C
C                  Derivatives of outside pressure w.r.t. ETA's
C                  (at fixed T)
      DPODEP = NPOUT*T+DV_DEP
      DPODEN = NNOUT*T+DV_DEN
C
C                  Derivatives of alpha density w.r.t. ETA's
C                  (at fixed T)
      DNADEP = ALFDNS*(2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP)/T
      DNADEN = ALFDNS*(2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN)/T
C
C                  Derivatives of alpha pressure w.r.t. ETA's
C                  (at fixed T)
      DPADEP = T*DNADEP
      DPADEN = T*DNADEN
C
C                  Derivatives of particle densities w.r.t. T
C                  (at fixed ETA's)
      DNPODT = 1.5*NPOUT/T
      DNNODT = 1.5*NNOUT/T
C
C                  Derivatives of exterior chem. pot. w.r.t. T
C                  (at fixed ETA's)
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT
C
C                  Derivative of pressure potential w.r.t. T
C                  (at fixed ETA's)
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT
C
C                  Derivative of outside pressure w.r.t. T
C                  (at fixed ETA's)
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT
C
C                  Derivative of alpha chem. pot. w.r.t. T
C                  (at fixed ETA's)
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT
C
C                  Derivative of alpha particle density w.r.t. T
C                  (at fixed ETA's)
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T
C
C                  Derivative of alpha particle pressure w.r.t. T
C                  (at fixed ETA's)
      DPADT = ALFDNS+T*DNADT
C
C
C                 ------------------------------------
C                 !      Derivatives of interior     !
C                 !      quantities                  !
C                 !      (w.r.t. Temp. and density)  !
C                 !                                  !
C                 ------------------------------------
C
C
C                   Derivatives of kinetic energy densities w.r.t. T
C                   (holding the number densities (X & NSUBI) fixed)
      DTPIDT =2.5*TAU_PI/T-2.25*X*NSUBI*GPI/UQ
      DTNIDT =2.5*TAU_NI/T-2.25*(1.0-X)*NSUBI*GNI/UQ
C
C                   Derivatives of pressures w.r.t. T
C                   (holding the number densities (X & NSUBI) fixed)
      DPIDT = OVR23*UQ*(DTPIDT+DTNIDT)
C
C                   Derivatives of interior chem. pot. w.r.t. T
C                   (holding the number densities (X & NSUBI) fixed)
      DMPIDT = ETA_PI-1.5*GPI
      DMNIDT = ETA_NI-1.5*GNI
C
C
C                  Derivatives of inside potentials w.r.t.
C                  interior proton and neutron densities
C                  (at fixed T)
c20      DVPIDP = 2.0*AA+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0))
      DVPIDP = DPVPDP(X*NSUBI,(1.0-X)*NSUBI)
c20      DVPIDN = 2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0))
      DVPIDN = DPVPDN(X*NSUBI,(1.0-X)*NSUBI)
c20      DVNIDP = DVPIDN
      DVNIDP = DPVNDP(X*NSUBI,(1.0-X)*NSUBI)
c20      DVNIDN = DVPIDP
      DVNIDN = DPVNDN(X*NSUBI,(1.0-X)*NSUBI)
C
C
C                   Derivatives of interior chemical potentials
C                   w.r.t. interior neutron and proton densities
C                  (at fixed T)
      DMPIDP = T*GPI/(X*NSUBI)+DVPIDP
      DMPIDN = DVPIDN
      DMNIDP = DVNIDP
      DMNIDN = T*GNI/((1.0-X)*NSUBI)+DVNIDN
C
C                   Derivatives of interior pressure
C                   w.r.t. interior neutron and proton densities
C                  (at fixed T)
      DPIDP = X*NSUBI*DMPIDP+(1.0-X)*NSUBI*DMNIDP
      DPIDN = X*NSUBI*DMPIDN+(1.0-X)*NSUBI*DMNIDN
C
C
C
C
C                 ------------------------------------
C                 !      Derivatives of "B" terms    !
C                 !      from the chemical and       !
C                 !      pressure equilibrium        !
C                 !      equations                   !
C                 !                                  !
C                 !      (w.r.t. Temperature )       !
C                 !                                  !
C                 ------------------------------------
C
C
C             Derivative of term from pressure equilibrium eqn.
C
      DB1DT = OVR23*ZETA*(SCRDUP-OVR23*SCRD)*HPRIM/H+
     1    ZETA*(SCRDPT-OVR23*SCRDT)-
     2    TRSCAL*U_NUC*NSUBI*(HPRIM*MUSUBT+H*DMUTDT)/AZERO
C
C
C             Derivative of term from proton equilibrium eqn.
C
      TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X
      TMP5 = DHDTDX/H-DHDX*HPRIM/H**2+
     1 1.5*SCRDXT/SCRDU-1.5*SCRDUX*SCRDUT/SCRDU**2
C
      DB2DT = OVR49*(ZETA*SCRD*HPRIM/(H*NSUBI))*TMP4+
     1    OVR23*ZETA*SCRDT*TMP4/NSUBI+
     2    OVR23*(ZETA*SCRD/NSUBI)*(X-1.0)*TMP5-
     3    TRSCAL*EXCLU*(DMUTDT*(H+DHDX*(1.0-X))+MUSUBT*
     4    (HPRIM+DHDTDX*(1.0-X))-DHDX*(1.0-X)-T*DHDX*(1.0-X))/AZERO
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C             Derivative of term from neutron equilibrium eqn.
C
      TMP4 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU
      TMP5 = DHDTDX/H-DHDX*HPRIM/H**2+
     1 1.5*SCRDXT/SCRDU-1.5*SCRDUX*SCRDUT/SCRDU**2
      DB3DT = OVR49*(ZETA*SCRD*HPRIM/(H*NSUBI))*X*TMP4+
     1        OVR23*(ZETA*SCRDT/NSUBI)*X*TMP4+
     2        OVR23*(ZETA*SCRD/NSUBI)*X*TMP5-
     3        TRSCAL*EXCLU*(HPRIM*MUSUBT+H*DMUTDT-X*DHDTDX*(MUSUBT-T)-
     4        X*DHDX*(DMUTDT-1.0))/AZERO
C
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
C                 ------------------------------------
C                 !      Derivatives of constraint   !
C                 !      and equilibrium equations   !
C                 !      with respect to the five    !
C                 !      compositional variables     !
C                 !      (U,x,n_i,eta_po,eta_no)     !
C                 !      and the three independent   !
C                 !      variables                   !
C                 !      (Baryon density, T, and Ye) !
C                 !                                  !
C                 ------------------------------------
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 1 (Baryon conservation)
C
      DFDOM(1,1) = NOUT*EXALFA+4.0*ALFDNS-NSUBI
C
      DFDOM(1,2) = 0.0
C
      DFDOM(1,3) = -U_NUC
C
      DFDOM(1,4) = -EXCLU*EXALFA*NPOUT/GPO+
     1             V_ALFA*DNADEP*EXCLU*NOUT-4.0*EXCLU*DNADEP
C
      DFDOM(1,5) = -EXCLU*EXALFA*NNOUT/GNO+
     1             V_ALFA*DNADEN*EXCLU*NOUT-4.0*EXCLU*DNADEN
C
C
C
      DFDL_1(1) = -1.0
C
      DFDL_2(1) = EXCLU*EXALFA*(DNPODT+DNNODT)-EXCLU*V_ALFA*NOUT*DNADT+
     1     4.0*EXCLU*DNADT
C
      DFDL_3(1) = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 2 (Charge conservation)
C
      DFDOM(2,1) = EXALFA*NPOUT+2.0*ALFDNS-X*NSUBI
C
      DFDOM(2,2) = -U_NUC*NSUBI
C
      DFDOM(2,3) = -X*U_NUC
C
      DFDOM(2,4) = -EXCLU*EXALFA*NPOUT/GPO+
     1     V_ALFA*EXCLU*NPOUT*DNADEP-2.0*EXCLU*DNADEP
C
      DFDOM(2,5) = V_ALFA*EXCLU*NPOUT*DNADEN-2.0*EXCLU*DNADEN
C
C
C
      DFDL_1(2) = -1.0*YE
C
      DFDL_2(2) = EXCLU*EXALFA*DNPODT-V_ALFA*EXCLU*NPOUT*DNADT+
     1     2.0*EXCLU*DNADT
C
      DFDL_3(2) = -1.0*BRYDNS
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 3 (Proton chemical equilibrium)
C
      DFDOM(3,1) = -DB2DU
C
      DFDOM(3,2) = NSUBI*(DMPIDP-DMPIDN)-DB2DX
C
      DFDOM(3,3) = (1.0-X)*DMPIDN+X*DMPIDP-DB2DNI
C
      DFDOM(3,4) = -DMPDEP
C
      DFDOM(3,5) = -DMPDEN
C
      DFDL_1(3) = 0.0
      DFDL_2(3) = -1.0*(DMPIDT-DMPODT-DB2DT)
      DFDL_3(3) = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 4 (Neutron chemical equilibrium)
C
      DFDOM(4,1) = -DB3DU
C
      DFDOM(4,2) = NSUBI*(DMNIDP-DMNIDN)-DB3DX
C
      DFDOM(4,3) = (1.0-X)*DMNIDN+X*DMNIDP-DB3DNI
C
      DFDOM(4,4) = -DMNDEP
C
      DFDOM(4,5) = -DMNDEN
C
      DFDL_1(4) = 0.0
      DFDL_2(4) = -1.0*(DMNIDT-DMNODT-DB3DT)
      DFDL_3(4) = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                Equation 5 (Pressure equilibrium)
C
      DFDOM(5,1) = -DB1DU
C
      DFDOM(5,2) = NSUBI*(DPIDP-DPIDN)-DB1DX
C
      DFDOM(5,3) = (1.0-X)*DPIDN+X*DPIDP-DB1DNI
C     ncomp = dfdom(5,3)
C
      DFDOM(5,4) = -DPODEP-DPADEP
C
      DFDOM(5,5) = -DPODEN-DPADEN
C
      DFDL_1(5) = 0.0
      DFDL_2(5) = -1.0*(DPIDT-DPODT-DPADT-DB1DT)
      DFDL_3(5) = 0.0
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
c      write(*,*) ' '
cc      write(*,*) ' '
cc      write(*,7010) db1dx,db2dx,db3dx
cc      write(*,7010) db1du,db2du,db3du
cc      write(*,7010) db1dni,db2dni,db3dni
cc      write(*,7010) db1dt,db2dt,db3dt
 7010 format(3(1x,g13.6))
c      write(*,7000) x,u_nuc,nsubi,eta_po,eta_no
c      write(*,*) 'eta_i ',eta_pi,eta_ni
c      write(*,*) ' as ',a1,a2,a3
c      write(*,*) ' '
c      write(*,7000) (dfdom(1,i),i=1,5,1)
c      write(*,7000) (dfdom(2,i),i=1,5,1)
c      write(*,7000) (dfdom(3,i),i=1,5,1)
c      write(*,7000) (dfdom(4,i),i=1,5,1)
c      write(*,7000) (dfdom(5,i),i=1,5,1)
c      write(*,*) ' '
cc      write(*,*) ' dna: ',dnadpo,dnadno
cc      write(*,*) ' dt: ',dmpidt,dmpodt
c      write(*,7000) (dfdl_1(i),i=1,5,1)
 7000 format(5(1x,g13.6))
c
c      pause
C                    Invert the DFDOM matrix
C
      CALL MATINV(DFDOM,DFDOMI,5)
C  IMSL subroutine call to invert the matrix
CCC      CALL DLINRG(5,DFDOM,5,DFDOMI,5)
C
cc      call matmul(dfdom,dfdomi,a_tmp,5,5,5)
c
cc      write(*,*) ' '
cc      write(*,7000) (a_tmp(1,i),i=1,5,1)
cc      write(*,7000) (a_tmp(2,i),i=1,5,1)
cc      write(*,7000) (a_tmp(3,i),i=1,5,1)
cc      write(*,7000) (a_tmp(4,i),i=1,5,1)
cc      write(*,7000) (a_tmp(5,i),i=1,5,1)
c
cc      write(*,7000) (dfdomi(1,i),i=1,5,1)
cc      write(*,7000) (dfdomi(2,i),i=1,5,1)
cc      write(*,7000) (dfdomi(3,i),i=1,5,1)
cc      write(*,7000) (dfdomi(4,i),i=1,5,1)
cc      write(*,7000) (dfdomi(5,i),i=1,5,1)
cc      write(*,*) ' >>>>>>>>>>>>>>>    dfdl_2 <<<<<<<<<<<<<<<<< '
cc      write(*,7000) (dfdl_2(i),i=1,5,1)
cc      write(*,*) ' >>>>>>>>>>>>>>>    dfdl_2 <<<<<<<<<<<<<<<<< '
c
c      DO 800 LLI=1,5,1
c        R_CHECK(LLI) = 0.0
c        DO 801 KKI=1,5,1
c          R_CHECK(LLI) = R_CHECK(LLI)+DFDOM(LLI,KKI)*RESULT(KKI)
c 801    CONTINUE
c        r_check(lli) = r_check(lli)-dfdl_1(lli)
c 800  CONTINUE
c      write(*,*) ' >>>>>>>>>>>>>>>    R check <<<<<<<<<<<<<<<<< '
c      write(*,7000) (r_check(i),i=1,5,1)
c
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                    Multiply the DFDL_1 vector by the DFDOMI matrix
C                    to get the density derivatives
C
      CALL MV_MUL(DFDOMI,DFDL_1,RESULT,5)
C
      DU_DN = RESULT(1)
      DX_DN = RESULT(2)
      DNI_DN = RESULT(3)
      DEP_DN = RESULT(4)
      DEN_DN = RESULT(5)
C
C
C                    Multiply the DFDL_2 vector by the DFDOMI matrix
C                    to get the Temperature derivatives
C
      CALL MV_MUL(DFDOMI,DFDL_2,RESULT,5)
C
      DU_DT = RESULT(1)
      DX_DT = RESULT(2)
      DNI_DT = RESULT(3)
      DEP_DT = RESULT(4)
      DEN_DT = RESULT(5)
C
C                    Multiply the DFDL_3 vector by the DFDOMI matrix
C                    to get the Ye derivatives
C
      CALL MV_MUL(DFDOMI,DFDL_3,RESULT,5)
C
      DU_DY = RESULT(1)
      DX_DY = RESULT(2)
      DNI_DY = RESULT(3)
      DEP_DY = RESULT(4)
      DEN_DY = RESULT(5)
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 ------------------------------------
C                 !      Derivatives of finite size  !
C                 !      terms in the internal       !
C                 !      energy and entropy          !
C                 !      densities w.r.t. to U,X,n_i !
C                 !      and T.  These are used in   !
C                 !      calculating the derivatives !
C                 !      w.r.t. the independant vars !
C                 !      (Baryon density, T, and Ye) !
C                 !                                  !
C                 ------------------------------------
C
C                        Free energy Surface & Coulomb terms
C                                  (Densities)
C
      F_SC = ZETA*SCRDU
C
      DFSCDU = ZETA*SCRDUP
C
      DFSCDX = ZETA*SCRDUX+SCRDU*DZDX
C
      DFSCDN = SCRDU*DZDNI
C
      DFSCDT = ZETA*SCRDUT+SCRDU*DZDT
C
C
C                        Free energy translational terms
C                                  (Densities)
      FTR = U_NUC*EXCLU*NSUBI*FTRANS
C
      DFTRDT = FTR*(HPRIM/H+1.0/T)-
     1    1.5*TRSCAL*U_NUC*EXCLU*NSUBI*H/AZERO
C
      DFTRDX = FTR*DHDX/H
C
      DFTRDU = FTR/U_NUC-FTR/EXCLU+
     1    TRSCAL*NSUBI*H*(1.0-2.0*U_NUC)/AZERO
C
      DFTRDN = FTR/NSUBI+TRSCAL*U_NUC*EXCLU*H*T/AZERO
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C                        Internal energy Surface & Coulomb terms
C                                  (Densities)
C
      TMP4 = 1.0-T*SCRDUT/SCRDU-OVR23*T*HPRIM/H
C
      E_SC = F_SC*TMP4
C
      DESCDU = DFSCDU*TMP4+
     1    F_SC*(T*SCRDUT*SCRDUP/SCRDU**2-T*SCRDPT/SCRDU)
C
      DESCDX = DFSCDX*TMP4+
     1    F_SC*(T*SCRDUT*SCRDUX/SCRDU**2-T*SCRDXT/SCRDU+
     2    OVR23*T*HPRIM*DHDX/H**2-OVR23*T*DHDTDX/H)
C
      DESCDN = DFSCDN*TMP4
C
      DESCDT = DFSCDT*TMP4+F_SC*
     1   (T*(SCRDUT**2)/SCRDU**2-SCRDUT/SCRDU-T*SCRDTT/SCRDU+
     2    OVR23*T*(HPRIM**2)/H**2-OVR23*HPRIM/H-OVR23*T*HPPRIM/H)
C
C                        Internal energy translational terms
C                                  (Densities)
C
      TMP4 = 1.5*H*T/AZERO-T*HPRIM*(MUSUBT-T)/AZERO
C
      E_TR = TRSCAL*EXCLU*BRYDNS*XH*TMP4
C
      DETRDU = TRSCAL*(NSUBI*(1.0-2.0*U_NUC)*TMP4-
     1    NSUBI*(T**2)*HPRIM*(1.0-2.0*U_NUC)/AZERO)
C
      DETRDX = TRSCAL*BRYDNS*XH*EXCLU*
     1    (1.5*T*DHDX/AZERO-T*(MUSUBT-T)*DHDTDX/AZERO)
C
      DETRDN = TRSCAL*(U_NUC*EXCLU*TMP4-
     1    BRYDNS*XH*EXCLU*(T**2)*HPRIM/(NSUBI*AZERO))
C
      DETRDT = TRSCAL*BRYDNS*XH*EXCLU*
     1    (1.5*(H+T*HPRIM)/AZERO-(HPRIM+T*HPPRIM)*(MUSUBT-T)/AZERO-
     2    T*HPRIM*(MUSUBT/T-2.5)/AZERO )
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C                        Entropy Surface & Coulomb terms
C                                  (Densities)
C
      S_SC = (E_SC-F_SC)/T
C
      DSSCDU = (DESCDU-DFSCDU)/T
C
      DSSCDX = (DESCDX-DFSCDX)/T
C
      DSSCDN = (DESCDN-DFSCDN)/T
C
      DSSCDT = (DESCDT-DFSCDT)/T-(E_SC-F_SC)/T**2
C
C                        Entropy translational terms
C                                  (Densities)
C
      TMP4 = MUSUBT*(HPRIM+H/T)/AZERO-(T*HPRIM+2.5*H)/AZERO
C
      S_TR = -TRSCAL*BRYDNS*XH*EXCLU*TMP4
C
      DSTRDU = -TRSCAL*(NSUBI*(1.0-2.0*U_NUC)*TMP4+
     1    NSUBI*T*(1.0-2.0*U_NUC)*(HPRIM+H/T)/AZERO)
C
      DSTRDX = -TRSCAL*BRYDNS*XH*EXCLU*
     1    (MUSUBT*(DHDTDX+DHDX/T)/AZERO-
     2    (T*DHDTDX+2.5*DHDX)/AZERO)
C
      DSTRDN = -TRSCAL*
     1    (U_NUC*EXCLU*TMP4+U_NUC*EXCLU*T*(HPRIM+H/T)/AZERO)
C
      DSTRDT = -(BRYDNS*XH*EXCLU*((MUSUBT/T-1.5)*(HPRIM+H/T)/AZERO+
     1    MUSUBT*(HPPRIM+HPRIM/T-H/T**2)/AZERO-
     2    (3.5*HPRIM+T*HPPRIM)/AZERO ))*TRSCAL
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of interior bulk !
C                 !      terms in the internal        !
C                 !      energy and entropy           !
C                 !      densities w.r.t. to U,X,n_i  !
C                 !      and T.  These are used in    !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
C
      S_NUC =(OVR53*UQ/T)*(TAU_NI+TAU_PI)-
     1    NSUBI*((1.0-X)*ETA_NI+X*ETA_PI)
C
c20      E_NUC = UQ*(TAU_PI+TAU_NI)+(NSUBI**2)*(AA+4.0*BB*X*(1.0-X))+
c20     1    CC*NSUBI**(1.0+DD)+X*NSUBI*DELTAM
      E_NUC = UQ*(TAU_PI+TAU_NI)+PV_E(X*NSUBI,(1.0-X)*NSUBI)
C
C
C                    Interior particle densties
      NPI = X*NSUBI
      NNI = (1.0-X)*NSUBI
C
      DTPIDT = 2.5*TAU_PI/T-2.25*NPI*GPI/UQ
      DTNIDT = 2.5*TAU_NI/T-2.25*NNI*GNI/UQ
C
C               Derivative of interior entropy density w.r.t. T
      DSIDT = UQ*(DTPIDT+DTNIDT)/T
C
C               Derivative of interior internal energy density w.r.t. T
      DEIDT = T*DSIDT
C
C
C
C
C                    Derivatives of eta's w.r.t. X and NSUBI
      DETPDX = GPI/X
      DETNDX = -GNI/(1.0-X)
      DETPDN = GPI/NSUBI
      DETNDN = GNI/NSUBI
C
C                    Derivatives of Tau's w.r.t. X and NSUBI
      DTPIDX = 1.5*T*NPI*DETPDX/UQ
      DTNIDX = 1.5*T*NNI*DETNDX/UQ
      DTPDNI = 1.5*T*NPI*DETPDN/UQ
      DTNDNI = 1.5*T*NNI*DETNDN/UQ
C
C
C
C           Derivative of interior entropy density w.r.t. X
      DSIDX = OVR53*UQ*(DTPIDX+DTNIDX)/T-NSUBI*(ETA_PI-ETA_NI)-
     1    NSUBI*((1.0-X)*DETNDX+X*DETPDX)
C
C           Derivative of interior internal energy density w.r.t. X
c20      DEIDX = UQ*(DTPIDX+DTNIDX)+
c20     1    (NSUBI**2)*4.0*BB*(1.0-2.0*X)+NSUBI*DELTAM
      DEIDX = UQ*(DTPIDX+DTNIDX)+DPVEDX(NSUBI,X)
C
C
C           Derivative of interior entropy density w.r.t. NSUBI
      DSIDN = OVR53*UQ*(DTPDNI+DTNDNI)/T-((1.0-X)*ETA_NI+X*ETA_PI)-
     1    NSUBI*((1.0-X)*DETNDN+X*DETPDN)
C
C
C           Derivative of interior internal energy density w.r.t. NSUBI
c20      DEIDN = UQ*(DTPDNI+DTNDNI)+2.0*NSUBI*(AA+4.0*BB*X*(1.0-X))+
c20     1    CC*(1.0+DD)*(NSUBI**DD)+X*DELTAM
      DEIDN = UQ*(DTPDNI+DTNDNI)+DPVEDN(NSUBI,X)
C
C
C
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of exterior bulk !
C                 !      nucleon internal energy &    !
C                 !      entropy densities and the    !
C                 !      chem. pot.  w.r.t. to eta_p, !
C                 !      ate_n & T. These are used in !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
      S_OUT =(OVR53*UQ/T)*(TAU_NO+TAU_PO)-NNOUT*ETA_NO-NPOUT*ETA_PO
C
c20      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*
c20     1((NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM)
      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*PV_E(NPOUT,NNOUT)
C
C                   Derivative of exterior entropy density w.r.t. T
      DSODT =  OVR53*UQ*(1.5*(TAU_PO+TAU_NO)/(T**2))-
     1     1.5*(NPOUT*ETA_PO+NNOUT*ETA_NO)/T
C
      DEODT = T*DSODT
C
C                    Derivatives of exterior particle densities w.r.t.
C                    Temperature (ETA's fixed)
      DNPODT = 1.5*NPOUT/T
      DNNODT = 1.5*NNOUT/T
C
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT
C
C
      DNPDEP = NPOUT/GPO
      DNNDEN = NNOUT/GNO
C
      DTPDEP = 1.5*T*NPOUT/UQ
      DTNDEN = 1.5*T*NNOUT/UQ
C
      DSODEP = (OVR53*UQ/T)*DTPDEP-NPOUT-ETA_PO*DNPDEP
      DSODEN = (OVR53*UQ/T)*DTNDEN-NNOUT-ETA_NO*DNNDEN
C
C
C                    Exterior particle potentials
c20      VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD )
      VNOUT = EIFLAG*PVN(NPOUT,NNOUT)
c20      VPOUT = EIFLAG*
c20     1    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*(1.0+DD)*NOUT**DD+DELTAM)
      VPOUT = EIFLAG*PVP(NPOUT,NNOUT)
C
C
      DEODEP = UQ*DTPDEP+VPOUT*DNPDEP
      DEODEN = UQ*DTNDEN+VNOUT*DNNDEN
C
      DMPDEP = T+DVPODP*NPOUT/GPO
      DMPDEN = DVPODN*NNOUT/GNO
      DMNDEP = DVNODP*NPOUT/GPO
      DMNDEN = T+DVNODN*NNOUT/GNO
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C                 -------------------------------------
C                 !      Derivatives of alpha         !
C                 !      particle internal energy &   !
C                 !      entropy densities and the    !
C                 !      chem. pot.  w.r.t. to eta_p, !
C                 !      ate_n & T. These are used in !
C                 !      calculating the derivatives  !
C                 !      w.r.t. the independant vars  !
C                 !      (Baryon density, T, and Ye)  !
C                 !                                   !
C                 -------------------------------------
C
C
C
C
      S_ALFA = ALFDNS*(2.5-MUALFA/T)
C
C
      E_ALFA = ALFDNS*(1.5*T-BALPHA)
C
C                  Derivative of pressure potential w.r.t. T
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT
C
C                  Derivative of outside pressure w.r.t. T
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT
C
C
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT
C
C                  Derivative of alpha particle density w.r.t. T
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T
C
C
      DSADT = DNADT*(2.5-MUALFA/T)-ALFDNS*DMUADT/T+ALFDNS*MUALFA/T**2
C
      DEADT = DNADT*(1.5*T-BALPHA)+1.5*ALFDNS
C
C
      DV_DEP = DV_DPO*NPOUT/GPO
      DV_DEN = DV_DNO*NNOUT/GNO
C
      DPODEP = OVR23*UQ*DTPDEP+DV_DEP
      DPODEN = OVR23*UQ*DTNDEN+DV_DEN
C
      DMADEP = 2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP
      DMADEN = 2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN
C
      DNADEP = ALFDNS*DMADEP/T
      DNADEN = ALFDNS*DMADEN/T
C
      DSADEP = DNADEP*(2.5-MUALFA/T)-ALFDNS*DMADEP/T
      DSADEN = DNADEN*(2.5-MUALFA/T)-ALFDNS*DMADEN/T
C
      DEADEP = DNADEP*(1.5*T-BALPHA)
      DEADEN = DNADEN*(1.5*T-BALPHA)
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      S_DENS = U_NUC*S_NUC+EXCLU*EXALFA*S_OUT+EXCLU*S_ALFA+S_SC+S_TR
C
      E_DENS = U_NUC*E_NUC+EXCLU*EXALFA*E_OUT+EXCLU*E_ALFA+E_SC+E_TR
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !      Temperature Derivatives     !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
      DNA_DT = DNADT+DNADEP*DEP_DT+DNADEN*DEN_DT
C
C
      DBSDT = (DU_DT*S_NUC-
     1    DU_DT*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DT*S_OUT
     2    -DU_DT*S_ALFA+
     3    U_NUC*(DSIDT+DSIDX*DX_DT+DSIDN*DNI_DT)+
     4    EXCLU*EXALFA*(DSODT+DSODEP*DEP_DT+DSODEN*DEN_DT)+
     5    EXCLU*(DSADT+DSADEP*DEP_DT+DSADEN*DEN_DT)+
     6    DSSCDT+DSSCDU*DU_DT+DSSCDX*DX_DT+DSSCDN*DNI_DT+
     7    DSTRDT+DSTRDU*DU_DT+DSTRDX*DX_DT+DSTRDN*DNI_DT)/BRYDNS
C
C
C~~~~~~~~~~~~~~~~~~
C
      DBUDT = T*DBSDT
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDT = DBUDT-S_DENS/BRYDNS-T*DBSDT
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDT = YE*(DMPODT+DMPDEP*DEP_DT+DMPDEN*DEN_DT)+
     1    (1.0-YE)*(DMNODT+DMNDEP*DEP_DT+DMNDEN*DEN_DT)
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDT = BRYDNS*(DBMUDT-DBFDT)
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !       Density Derivatives        !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
C
      DNA_DN = DNADEP*DEP_DN+DNADEN*DEN_DN
C
C
      DBSDN = (DU_DN*S_NUC-
     1    DU_DN*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DN*S_OUT-DU_DN*S_ALFA+
     2    U_NUC*(DSIDX*DX_DN+DSIDN*DNI_DN)+
     3    EXCLU*EXALFA*(DSODEP*DEP_DN+DSODEN*DEN_DN)+
     4    EXCLU*(DSADEP*DEP_DN+DSADEN*DEN_DN)+
     5    DSSCDU*DU_DN+DSSCDX*DX_DN+DSSCDN*DNI_DN+
     6    DSTRDU*DU_DN+DSTRDX*DX_DN+DSTRDN*DNI_DN)/BRYDNS-
     7    S_DENS/BRYDNS**2
C
C
C
C~~~~~~~~~~~~~~~~~~
C
      DBUDN = (DU_DN*E_NUC-
     1    DU_DN*EXALFA*E_OUT-EXCLU*V_ALFA*DNA_DN*E_OUT-DU_DN*E_ALFA+
     2    U_NUC*(DEIDX*DX_DN+DEIDN*DNI_DN)+
     3    EXCLU*EXALFA*(DEODEP*DEP_DN+DEODEN*DEN_DN)+
     4    EXCLU*(DEADEP*DEP_DN+DEADEN*DEN_DN)+
     5    DESCDU*DU_DN+DESCDX*DX_DN+DESCDN*DNI_DN+
     6    DETRDU*DU_DN+DETRDX*DX_DN+DETRDN*DNI_DN)/BRYDNS-
     7    E_DENS/BRYDNS**2
C
C
C
C
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDN = DBUDN-T*DBSDN
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDN = YE*(DMPDEP*DEP_DN+DMPDEN*DEN_DN)+
     1    (1.0-YE)*(DMNDEP*DEP_DN+DMNDEN*DEN_DN)
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDN = BRYDNS*(DBMUDN-DBFDN)+MUBARY-BFTOT
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C
C                 ------------------------------------
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 !         Ye Derivatives           !
C                 !                                  !
C                 !                                  !
C                 !                                  !
C                 ------------------------------------
C
C
C
C
      DNA_DY = DNADEP*DEP_DY+DNADEN*DEN_DY
C
C
      DBSDY = (DU_DY*S_NUC-
     1    DU_DY*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DY*S_OUT-DU_DY*S_ALFA+
     2    U_NUC*(DSIDX*DX_DY+DSIDN*DNI_DY)+
     3    EXCLU*EXALFA*(DSODEP*DEP_DY+DSODEN*DEN_DY)+
     4    EXCLU*(DSADEP*DEP_DY+DSADEN*DEN_DY)+
     5    DSSCDU*DU_DY+DSSCDX*DX_DY+DSSCDN*DNI_DY+
     6    DSTRDU*DU_DY+DSTRDX*DX_DY+DSTRDN*DNI_DY)/BRYDNS
C
C
C~~~~~~~~~~~~~~~~~~
C
      DBUDY = (DU_DY*E_NUC-
     1    DU_DY*EXALFA*E_OUT-EXCLU*V_ALFA*DNA_DY*E_OUT-DU_DY*E_ALFA+
     2    U_NUC*(DEIDX*DX_DY+DEIDN*DNI_DY)+
     3    EXCLU*EXALFA*(DEODEP*DEP_DY+DEODEN*DEN_DY)+
     4    EXCLU*(DEADEP*DEP_DY+DEADEN*DEN_DY)+
     5    DESCDU*DU_DY+DESCDX*DX_DY+DESCDN*DNI_DY+
     6    DETRDU*DU_DY+DETRDX*DX_DY+DETRDN*DNI_DY)/BRYDNS
C
C
C
C~~~~~~~~~~~~~~~~~~
C
C
      DBFDY = DBUDY-T*DBSDY
C
C~~~~~~~~~~~~~~~~~~
C
      DBMUDY = YE*(DMPDEP*DEP_DY+DMPDEN*DEN_DY)+MUPROT+
     1    (1.0-YE)*(DMNDEP*DEP_DY+DMNDEN*DEN_DY)-MUN
C
C~~~~~~~~~~~~~~~~~~
C
      DBPDY = BRYDNS*(DBMUDY-DBFDY)
C
C
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C-----------------------------------------------------------------------
C                End of derivatives of thermodynamic variables
C-----------------------------------------------------------------------
C
C                  Total derivatives
C                  (Baryons+Electrons+Photons)
C
      DUDT = DBUDT+DEUDT+DPUDT
      DUDN = DBUDN+DEUDN+DPUDN
      DUDY = DBUDY+DEUDY+DPUDY
C
C
      DSDT = DBSDT+DESDT+DPSDT
      DSDN = DBSDN+DESDN+DPSDN
      DSDY = DBSDY+DESDY+DPSDY
C
C
      DPDT = DBPDT+DEPDT+DPPDT
      DPDN = DBPDN+DEPDN+DPPDN
      DPDY = DBPDY+DEPDY+DPPDY
C
C
      DMUDT = DBMUDT+YE*DEMUDT
      DMUDN = DBMUDN+YE*DEMUDN
      DMUDY = DBMUDY+YE*DEMUDY
C
C                Calculate the adiabatic index
      GAM_S = BRYDNS*DPDN/PTOT+T*(DPDT**2)/(BRYDNS*PTOT*DUDT)
C
C
C                Set the value of XPREV to X for use the next
C                time through
C
      XPREV = X
C
C                Save the value of the proton density to be used
C                by the "no nuclei" scheme on the next call
      P_PREV = NPOUT
C
C
C                Return the three internal compositional variables
      INPVAR(2) = NSUBI
      INPVAR(3) = ETA_PO
      INPVAR(4) = ETA_NO
C
  999 CONTINUE
C
C
      INPVARP(1)     = INPVAR(1)
      INPVARP(2)     = INPVAR(2)
      INPVARP(3)     = INPVAR(3)
      INPVARP(4)     = INPVAR(4)
      YEP            = YE
      BRYDNSP        = BRYDNS
      XPREVP         = XPREV
      P_PREVP        = P_PREV
      SSFLAGP        = SSFLAG
C
C                Rejoice for this routine is finished!!!!!!!
      RETURN
C
C
      END
Cnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnucnu

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         RESET.FOR
C
C***********************************************************************
C
C    MODULE:       RESET
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         12/21/90
C
C                  Please report any problems to me at:
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or
C                            fswesty@sbast3.sunysb.edu
C
C
C    CALL LINE:    CALL RESET(INPVAR,YE,BRYDNS,OUTVAR)
C
C
C    INPUTS:       INPVAR = TEMP, NSUBI, ETA_PO, ETA_NO
C                  YE = ELECTRON FRACTION
C                  BRYDNS = BARYON NUMBER DENSITY
C
C
C
C    OUTPUTS:      OUTVAR = ARRAY OF LENGTH 4 CONTAINING RESET VALUES
C                  FOR THE INITIAL GUESSES
C
C
C
C
C    INCLUDE FILES: NONE
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE RESET(INPVAR,YE,BRYDNS,OUTVAR)
C
C
      IMPLICIT NONE
C
C
C                      Subroutine parameters
C
      DOUBLE PRECISION INPVAR(4), OUTVAR(4), YE, BRYDNS
C
C
C                      Local variables
C
      DOUBLE PRECISION ZPG, ZNG, ETA_PG, ETA_NG, PI, UQ, MQ, T, EFRAC
C
C                      Functions
C
      DOUBLE PRECISION FINV12
C
C-----------------------------------------------------------------------
C
      T = INPVAR(1)
C
C
      PI = 3.1415927
      UQ = 20.721
C
      MQ = (T/UQ)**1.5
C
C
      EFRAC = 0.5*YE
C
      ZNG = 2.0*(PI**2)*BRYDNS*(1.0-EFRAC)/MQ
C
      ZPG = 2.0*(PI**2)*BRYDNS*EFRAC/MQ
C
      ETA_NG = FINV12(ZNG)
C
      ETA_PG = FINV12(ZPG)
C
      OUTVAR(1) = INPVAR(1)
      OUTVAR(2) = INPVAR(2)
      OUTVAR(3) = ETA_PG
      OUTVAR(4) = ETA_NG
C
C
C-----------------------------------------------------------------------
C
 999  RETURN
C
C
      END

!///////////////////////////////////////////////////////////////////////

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         SILEOS.F
C
C***********************************************************************
C
C    MODULE:       SILEOS
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: swesty@nuclear.physics.sunysb.edu
C                            -or-
C                            fswesty@sbast3.ess.sunysb.edu
C
C    VERSION:      1.0 1/20/93
C                  2.0 1/28/93
C
C
C    PURPOSE:      Calculate the EOS assuming a Boltzmann gas of Si^28
C                  plus electrons
C
C
C
C    CALL LINE:    CALL SILEOS(MODE,I_T,YE,BRYDNS,EB)
C
C    INPUTS:       MODE = Input variable mode (INTEGER)
C                        1--> Temp; 2--> Internal energy per baryon
C                  I_T = Temperature (in MeV) (DOUBLE PRECISION)
C                        (Used if MODE = 2)
C                  YE = Electron fraction (DOUBLE PRECISION)
C                  BRYDNS = Baryon density (in fm^-3) (DOUBLE PRECISION)
C                  EB = Internal energy per baryon (DOUBLE PRECISION)
C                       in units of MeV. (Used if MODE = 2)
C
C    OUTPUTS:      NONE
C
C
C    INCLUDE FILES:  EOS_M4A.INC, EL_EOS.INC
C
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SILEOS(MODE,TP,YEP,BRYDNSP,EB)
C
      USE el_eos_module
      USE eos_m4c_module
C
      IMPLICIT NONE
C
      INTEGER MODE
      DOUBLE PRECISION EB, DEL_T
C
C
C                      Baryon number density
      DOUBLE PRECISION N_SI, MU_SI
      DOUBLE PRECISION TP,YEP,BRYDNSP
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C                      Temperature as input
C23456789012345678901234567890123456789012345678901234567890123456789012
      T              = TP
      YE             = YEP
      BRYDNS         = BRYDNSP
C
      IF(MODE.EQ.1) THEN
C                      Number density of silicon nuclei
        N_SI = BRYDNS/28.0D0
C                      Atomic mass
        A = 28.0D0
C                      Proton fraction of Si^28
        X = 0.5D0
C                      Mass fractions
        XH = 1.0D0
        XALFA = 0.0D0
        XNUT = 0.0D0
        XPROT = 0.0D0
C
C                      Quantum concentration (Note the factor of 28
C                      is from the mass of silicon)
        NQ = 2.36D-4*((28.0D0*T)**1.5D0)
C
C                      Translational chemical potential of silicon
        MU_SI = T*DLOG(N_SI/NQ)
C
C                      Entropy per baryon of silicon
        BS = 2.5D0-MU_SI/T
C
C                      Silicon internal energy
        BU = 1.5D0*T-8.791
C
C                      Free energy per baryon
        BFTOT = BU-T*BS
C
C                      Silicon pressure
        BPRESS = N_SI*T
C
C                      Adiabatic index
        GAM_S = 1.667D0
C
C                      Difference in chemical potentials
        MUHAT = 0.0D0
C
C                      Call the electron EOS
        CALL EL_EOS(T,YE,BRYDNS)
C
C                      Total free energy per baryon
        FTOT = BFTOT+FSUBE+PF
C
C                      Total internal energy per baryon
        UTOT = BU+EU+PU
C
C                      Total entropy per baryon
        STOT = BS+ES+PS
C
C                      Total pressure
        PTOT = BPRESS+EPRESS+PPRESS
C
C
C                      Temperature derivatives
        DUDT = 1.5D0+DEUDT+DPUDT
        DPDT = N_SI+DEPDT+DPPDT
C                      Density derivatives
        DUDN = 0.0D0+DEUDN+DPUDN
        DPDN = T+DEPDN+DPPDN
C                      Ye derivatives
        DUDY = 0.0D0+DEUDY+DPUDY
        DPDY = 0.0D0+DEPDY+DPPDY
      ELSEIF(MODE.EQ.2) THEN
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C                      Internal energy per baryon as input
C23456789012345678901234567890123456789012345678901234567890123456789012
        DO 20 K=1,20,1
C                      Number density of silicon nuclei
          N_SI = BRYDNS/28.0D0
C                      Atomic mass
          A = 28.0D0
C                      Proton fraction of Si^28
          X = 0.5D0
C                      Mass fractions
          XH = 1.0D0
          XALFA = 0.0D0
          XNUT = 0.0D0
          XPROT = 0.0D0
C
C                      Quantum concentration (Note the factor of 28
C                      is from the mass of silicon)
          NQ = 2.36D-4*((28.0D0*T)**1.5D0)
C
C                      Translational chemical potential of silicon
          MU_SI = T*DLOG(N_SI/NQ)
C
C                      Entropy per baryon of silicon
          BS = 2.5D0-MU_SI/T
C
C                      Silicon internal energy
          BU = 1.5D0*T-8.791
C
C                      Free energy per baryon
          BFTOT = BU-T*BS
C
C                      Silicon pressure
          BPRESS = N_SI*T
C
C                      Adiabatic index
          GAM_S = 1.667D0
C
C                      Difference in chemical potentials
          MUHAT = 0.0D0
C
C                      Call the electron EOS
          CALL EL_EOS(T,YE,BRYDNS)
C
C                      Total free energy per baryon
          FTOT = BFTOT+FSUBE+PF
C
C                      Total internal energy per baryon
          UTOT = BU+EU+PU
C
C                      Total entropy per baryon
          STOT = BS+ES+PS
C
C                      Total pressure
          PTOT = BPRESS+EPRESS+PPRESS
C
C
C                      Temperature derivatives
          DUDT = 1.5D0+DEUDT+DPUDT
          DPDT = N_SI+DEPDT+DPPDT
C                      Density derivatives
          DUDN = 0.0D0+DEUDN+DPUDN
          DPDN = T+DEPDN+DPPDN
C                      Ye derivatives
          DUDY = 0.0D0+DEUDY+DPUDY
          DPDY = 0.0D0+DEPDY+DPPDY
C
          DEL_T = (EB-UTOT)/DUDT
C
          IF(DABS(DEL_T/T).LT.1.0D-6) THEN
            GOTO 30
          ELSE
            T = T+DEL_T
          ENDIF
 20     CONTINUE
C
 30     CONTINUE
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C                      Scream at the user & die!!!
C23456789012345678901234567890123456789012345678901234567890123456789012
      ELSE
        WRITE(*,*) ' SILEOS: Mode # ',MODE,' not supported! '
        STOP
C
      ENDIF
C
  999 CONTINUE
C
      TP             = T
      YEP            = YE
      BRYDNSP        = BRYDNS
C
C
      RETURN
      END

!///////////////////////////////////////////////////////////////////////

        SUBROUTINE SPLINE(X,Y,N,Y2)
c  Computes spline coefficients; Y(X) is input data; Y2 is output.
        implicit real*8(a-h,o-z)
        DIMENSION X(N),Y(N),Y2(N),U(500)
        Y2(1)=0.
        U(1)=0.
        DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
 11     U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     1 /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
        Y2(N)=0.
        DO 12 K=N-1,1,-1
 12     Y2(K)=Y2(K)*Y2(K+1)+U(K)
        RETURN
        END

!///////////////////////////////////////////////////////////////////////

        SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,KLO,KHI)
c     Computes spline fit of Y(X); YA(XA) is input data, Y2A are spline
c  coefficents, klo and khi are running indices which bracket X.
        implicit real*8(a-h,o-z)
        DIMENSION XA(N),YA(N),Y2A(N)
cc  Determine the bracketing indices
        IF((KHI.GT.N).OR.(KLO.LT.1)) THEN
          WRITE(*,*) ' SPLINT CRASHED '
          WRITE(*,*) ' N,KHI,KLO =',N,KHI,KLO
          STOP
        ENDIF
C
        IF((KHI.GT.N).OR.(KHI.LT.1)) KHI = N
        IF((KLO.GT.N).OR.(KLO.LT.1)) KLO = 1
C
        IF(KHI-KLO.GT.1)GOTO 1
c-----------------------------------------------------
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2
        KHI=MAX(KHI-1,1)
        KLO=MAX(KLO-1,1)
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2
        KHI=MIN(KHI+2,N)
        KLO=MIN(KLO+2,N)
        IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2
c-----------------------------------------------------
        KLO=1
        KHI=N
 1      IF(KHI-KLO.EQ.1) GOTO 2
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
        KHI=K
        ELSE
        KLO=K
        ENDIF
        GOTO 1
 2      H=XA(KHI)-XA(KLO)
        IF(H.EQ.0.) PAUSE 'BAD XA INPUT. '
cc  Compute spline fit.
        A=(XA(KHI)-X)/H
        B=(X-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+
     1 ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*H**2/6.
c       write(*,5)klo,khi,x,xa(klo),xa(khi),ya(klo),ya(khi)
c     > ,y2a(klo),y2a(khi),y
 5      format(2i3,8(1pe9.2))
        RETURN
        END

!///////////////////////////////////////////////////////////////////////

C
C
C***********************************************************************
C
C    MODULE:       VECCOP
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         4/10/90
C
C    PURPOSE:      Copy one vector on length N into another
C
C    CALL LINE:    CALL VECCOP(A,B,N)
C
C    INPUTS:       A = Array to be copied  (D)
C                  N = Number of rows in arrays (I)
C
C    OUTPUTS:      B = Array to be copied into (D)
C
C    CALLS :       None
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE VECCOP(A,B,N)
C
      IMPLICIT NONE
C
      INTEGER N
      DOUBLE PRECISION A(N), B(N)
C
C
C                 Local variables
C
      INTEGER I
C
      DO 10 I=1,N,1
        B(I) = A(I)
 10   CONTINUE
C
 999  RETURN
C
      END

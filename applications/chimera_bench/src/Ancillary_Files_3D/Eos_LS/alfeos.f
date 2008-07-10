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

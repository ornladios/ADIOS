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

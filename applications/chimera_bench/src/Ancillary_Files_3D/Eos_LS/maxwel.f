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

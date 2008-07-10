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

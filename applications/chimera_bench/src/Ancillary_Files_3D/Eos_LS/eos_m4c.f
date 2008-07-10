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

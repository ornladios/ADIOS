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

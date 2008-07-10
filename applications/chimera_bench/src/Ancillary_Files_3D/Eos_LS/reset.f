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

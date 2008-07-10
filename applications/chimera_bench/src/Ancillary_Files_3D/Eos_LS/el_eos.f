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

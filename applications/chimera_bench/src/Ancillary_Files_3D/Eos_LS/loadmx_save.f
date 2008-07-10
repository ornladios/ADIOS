C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         ALOADMX.FOR
C    MODULE:       LOADMX
C    TYPE:         LOADMX
C
C    PURPOSE:      LOAD THE LOOK-UP TABLE FOR THE MAXWELL CONSTRUCTION
C
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         7/16/90
C
C    CALL LINE:    CALL LOADMX
C
C    INPUTS:       N/A
C
C    OUTPUTS       N/A
C
C    SUBROUTINE CALLS: EOS_M4C
C
C    INCLUDE FILES:  EOS_M4C.INC, MAXWEL.INC
C
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE LOADMX()
C
C
      USE eos_m4c_module
      USE maxwel_module
C
      IMPLICIT NONE
C
      INTEGER NTMP, NYE, NYE2, NUM_BP
      INTEGER LUN1, LUN2, KK, KMIN
      PARAMETER(LUN1=54,LUN2=55)
C
C
      INTEGER FNML1, FNML2
      CHARACTER*60 FNAME1, FNAME2
      logical loaddata
C
      DOUBLE PRECISION N_SM, SYMM_M, COMP_M, BINDEM, SYM_SM, SIG_SM
      DOUBLE PRECISION N_SB, SYMM_B, COMP_B, BINDEB, SYM_SB, SIG_SB
      DOUBLE PRECISION MSCDD3, BSCDD3
C
      INCLUDE 'force.inc'
      data loaddata /.false./
C
C
C
c     CALL GETFNM(FNAME1,FNML1,'Enter ASCII Maxwell fname:',26)
C
c     CALL GETFNM(FNAME2,FNML2,'Enter ASCII boundary fname:',27)
          if ( loaddata ) return
      FNAME1        = '../../Eos_ls/max180.atb'
      FNML1         = 23
      FNAME2        = '../../Eos_ls/bd180.atb'
      FNML2         = 22
      loaddata      = .true.
C
C
C
C-----------------------------------------------------------------------
C        Read the file Maxwell construction data file
C-----------------------------------------------------------------------
C
C
C
      OPEN(UNIT=LUN1,FILE=FNAME1(1:FNML1),STATUS='OLD')
C
C
C
C
C
C
      READ(LUN1,*) N_SM, SYMM_M
      READ(LUN1,*) COMP_M,BINDEM
      READ(LUN1,*) SYM_SM, SIG_SM
C
C
C
C
      READ(LUN1,*) NTMP,NYE
      READ(LUN1,*) T_LOW,T_HI
      READ(LUN1,*) Y_LOW,Y_HI
C
C
C
      IF((NTMP.NE.NUMTMP).OR.(NYE.NE.NUMYE)) THEN
        WRITE(*,*) 'LOADMX:  MXWL TABLE IS INCOMPATIBLE WITH ARRAYS'
        STOP
      ENDIF
C
C
      DO 101 J=1,NUMYE,1
        DO 100 I=1,NUMTMP,3
          KMIN = MIN0(I+2,NUMTMP)
          READ(LUN1,*) (BRYLOW(KK,J),KK=I,KMIN,1)
 100    CONTINUE
 101  CONTINUE
C
C
      DO 103 J=1,NUMYE,1
        DO 102 I=1,NUMTMP,3
          KMIN = MIN0(I+2,NUMTMP)
          READ(LUN1,*) (BRYHI(KK,J),KK=I,KMIN,1)
 102    CONTINUE
 103  CONTINUE
C
C
C
      DO 104 I=1,NUMYE,3
        KMIN = MIN0(I+2,NUMYE)
        READ(LUN1,*) (T_H(KK),KK=I,KMIN,1)
 104  CONTINUE
C
C
      DO 105 I=1,NUMYE,3
        KMIN = MIN0(I+2,NUMYE)
        READ(LUN1,*) (D_H(KK),KK=I,KMIN,1)
 105  CONTINUE
C
      READ(LUN1,*) YCUT
      READ(LUN1,*) MSCDD3
C
C
      CLOSE(UNIT=LUN1,STATUS='KEEP')
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  MAXWELL CON. TABLE IS INITIALIZED>>'
      WRITE(*,*)
C
C
C
C-----------------------------------------------------------------------
C        Read the file Boundary data file
C-----------------------------------------------------------------------
C
C
C
      OPEN(UNIT=LUN2,FILE=FNAME2(1:FNML2),STATUS='OLD')
C
C
C
C
      READ(LUN2,*) N_SB,SYMM_B
      READ(LUN2,*) COMP_B,BINDEB
      READ(LUN2,*) SYM_SB,SIG_SB
C
C
C
C
C
C
      READ(LUN2,*) NUM_BP,NYE2
      READ(LUN2,*) LNL,LNH,LNC
      READ(LUN2,*) Y_LOW2,Y_HI2
C
C
      IF((NBPNTS.NE.NUM_BP).OR.(NYE2.NE.NUMYE)) THEN
        WRITE(*,*) 'LOADMX:  BNDY TABLE IS INCOMPATIBLE WITH ARRAYS'
        STOP
      ENDIF
C
      IF(ABS(LNL-LNLOW).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  LOWER END OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
C
      IF(ABS(LNH-LNHI).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  UPPER END OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
C
      IF(ABS(LNC-LNCUT).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  MID CUT OF PHASE BNDY IS INCONSIST.'
        STOP
      ENDIF
C
      IF(ABS(Y_LOW-Y_LOW2).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  LOWER YE LIMITS ARE INCONSIST.'
        STOP
      ENDIF
C
      IF(ABS(Y_HI-Y_HI2).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  UPPER YE LIMITS ARE INCONSIST.'
        STOP
      ENDIF
C
C
      DO 201 J=1,NUMYE,1
        DO 200 I=1,NBPNTS,3
          KMIN = MIN0(I+2,NBPNTS)
          READ(LUN2,*) (LBOUND(KK,J),KK=I,KMIN,1)
 200    CONTINUE
 201  CONTINUE
C
C
      DO 203 J=1,NUMYE,1
        DO 202 I=1,NBPNTS,3
          KMIN = MIN0(I+2,NBPNTS)
          READ(LUN2,*) (UBOUND(KK,J),KK=I,KMIN,1)
 202    CONTINUE
 203  CONTINUE
C
      READ(LUN2,*) BSCDD3
C
      IF(ABS(MSCDD3-BSCDD3).GT.1.0D-10) THEN
        WRITE(*,*) 'LOADMX:  SCRDD3 VALUES ARE INCONSIST.'
        STOP
      ENDIF
C
C
C
      CLOSE(UNIT=LUN2,STATUS='KEEP')
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  BOUNDARY TABLE IS INITIALIZED>>'
      WRITE(*,*)
C
C
C
C-----------------------------------------------------------------------
C                  All arrays are now loaded so return
C-----------------------------------------------------------------------
C
      SCRDD3 = BSCDD3
      N_S = N_SM
      NSUBS = N_SM
      SYMM = SYMM_M
      COMP = COMP_M
      BIND_E = BINDEM
      SYM_S = SYM_SM
      SIG_S = SIG_SM
C
c20      SKYRMC=(.3*((HBAR*C)**2)/MASSN)*(1.5*N_S*(PI**2))**OVR23
c20      DD = (COMP+2.0*SKYRMC)/(3.0*SKYRMC+9.0*BIND_E)
c20      BB = (SKYRMC*(2.0**OVR23-1.0)-SYMM)/N_S
c20      AA = (OVR23*SKYRMC-DD*(SKYRMC+BIND_E))/(N_S*(DD-1.0))-BB
c20      CC = (COMP+2.0*SKYRMC)/(9.0*DD*(DD-1.0)*N_S**DD)
      SKYRMC=(.3*((HBAR*C)**2)/MASSN)*(1.5*N_S*(PI**2))**OVR23
      DD = (COMP+2.0*SKYRMC+SCRDD3*(COMP-4.0*SKYRMC-18.0*BIND_E))/
     1    ((1.0-SCRDD3)*(3.0*SKYRMC+9.0*BIND_E))
      BB = (SKYRMC*(2.0**OVR23-1.0)-SYMM)/N_S
      CC = ((OVR3*SKYRMC+BIND_E)*(1.0+SCRDD3)**2)/((N_S**DD)*(DD-1.0))
      AA = ((OVR23*SKYRMC-DD*(SKYRMC+BIND_E)-SCRDD3*(OVR3*SKYRMC+BIND_E)
     1    )/(N_S*(DD-1.0)) )-BB
      DD3 = SCRDD3/(N_S**(DD-1.0))
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX:  SKYRME PARAMETERS FOR THIS RUN ARE:>>'
      WRITE(*,*) 'ABCD: ',AA,BB,CC,DD,SCRDD3
      WRITE(*,*) ' Satur. density, symmetry engy, & compression mod.:'
      WRITE(*,*) N_SM, SYMM_M, COMP_M
      WRITE(*,*) N_SB, SYMM_B, COMP_B
      WRITE(*,*) ' Binding engy, surf. symm. engy, & surface tension:'
      WRITE(*,*) BINDEM,SYM_SM,SIG_SM
      WRITE(*,*) BINDEB,SYM_SB,SIG_SB
C
      WRITE(*,*)
C
C
      CALL INITFERM()
C
      WRITE(*,*)
      WRITE(*,*) '<<LOADMX: FERMI INTEGRAL TABLES ARE INITIALIZED>>'
      WRITE(*,*)
C
C
 999  RETURN
C
      END

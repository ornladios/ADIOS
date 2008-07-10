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

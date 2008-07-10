C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    FILE:         EL_GL.FOR
C
C***********************************************************************
C
C    MODULE:       EL_GL
C    TYPE:         SUBROUTINE
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C
C    DATE:         2/12/91
C                  version 2.0 1/24/93 (includes full fermi integrals)
C
C                  BITNET:  SWESTY@SUNYSBNP or
C                  INTERNET: swesty@nuclear.physics.sunysb.edu or
C                            fswesty@sbast3.ess.sunysb.edu
C
C    PURPOSE:      The elctron and photon equation of state
C
C
C    CALL LINE:    CALL EL_GL(T,YE,BRYDNS)
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
      SUBROUTINE EL_GL(T,YE,BRYDNS)
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
      DOUBLE PRECISION MASS_E
      PARAMETER(MASS_E=0.51100D0)
C
C                           Pi and 1/3
      DOUBLE PRECISION PI, PI2, OVR3, MOVR3, OVR23
      PARAMETER(PI=3.1415927, PI2=9.8696044)
      PARAMETER(OVR3=0.33333333, MOVR3=-0.33333333, OVR23=0.66666667)
C
C                           2nd Fermi integral
      DOUBLE PRECISION F_2
C                           Positron degeneracy parameter
      DOUBLE PRECISION ELPETA
C                           Multiplicative factors in the G-L
C                           quadrature
      DOUBLE PRECISION TFAC, DTFDT, EFAC, PFAC
C                           Electron eta and beta
      DOUBLE PRECISION ETA_E, BETA, DBTDT, DETDT, DETDNE, D_ETA
C                           Number density and it's derivatives
      DOUBLE PRECISION NE_CHK, NEN, DNENDB, DNENDE
C                           Energy density & derivatives
      DOUBLE PRECISION E_ENG, DEDB, DEDE
C                           Pressure & derivatives
      DOUBLE PRECISION E_PR, DPDB, DPDE
C                           Electron & positron occupation numbers
      DOUBLE PRECISION FE, FP, DFEDB, DFEDE, DFPDB, DFPDE
C                           Integrand factors & derivatives
      DOUBLE PRECISION XN1, DXN1DB, XROOT
      DOUBLE PRECISION XN2, DXN2DB, XN3, DXN3DB
C                           Loop variables
      INTEGER I, J, MAXITR
      PARAMETER(MAXITR=30)
C
C                           Convergence tolerance for N-R iteration
      DOUBLE PRECISION EPSIL
      PARAMETER(EPSIL=1.0D-10)
C
      INTEGER NGAUSS
      PARAMETER(NGAUSS=16)
c      PARAMETER(NGAUSS=36)
      DOUBLE PRECISION WGHT(NGAUSS), X(NGAUSS)
C
      DATA X /.087649410479D00
     *,.46269632892D00,.11410577748D01,.21292836451D01,.34370866339D01
     *,.50780186145D01,.70703385350D01,.94383143364D01,.12214223369D02
     *,.15441527369D02,.19180156857D02,.23515905694D02,.28578729743D02
     *,.34583398702D02,.41940452648D02,.51701160340D02/
      DATA WGHT /.22503631486D00
     *,.52583605276D00,.83196139169D00,.11460992410D01,.14717513170D01
     *,.18131346874D01,.21755175197D01,.25657627502D01,.29932150864D01
     *,.34712344831D01,.40200440864D01,.46725166077D01,.54874206580D01
     *,.65853612333D01,.82763579844D01,.11824277552D02/
C
C23456789012345678901234567890123456789012345678901234567890123456789012
c      DATA X /
c     * 3.961320640860D-02 ,  2.088002856543D-01 ,  5.135107755352D-01 ,
c     * 9.543811530791D-01 ,  1.532249225394D+00 ,  2.248215817854D+00 ,
c     * 3.103661490498D+00 ,  4.100262545699D+00 ,  5.240010108755D+00 ,
c     * 6.525233330468D+00 ,  7.958627508999D+00 ,  9.543288035279D+00 ,
c     * 1.128275128666D+01 ,  1.318104390025D+01 ,  1.524274226656D+01 ,
c     * 1.747304463141D+01 ,  1.987785893308D+01 ,  2.246391051439D+01 ,
c     * 2.523887525730D+01 ,  2.821154567688D+01 ,  3.139204037421D+01 ,
c     * 3.479207144800D+01 ,  3.842529076480D+01 ,  4.230774567208D+01 ,
c     * 4.645849004290D+01 ,  5.090042150332D+01 ,  5.566145791468D+01 ,
c     * 6.077624068453D+01 ,  6.628869068025D+01 ,  7.225601475446D+01 ,
c     * 7.875533816278D+01 ,  8.589548142950D+01 ,  9.383992978484D+01 ,
c     * 1.028580101471D+02 ,  1.134651354897D+02 ,  1.269880151966D+02 /
C
c      DATA WGHT /
c     * 1.016694751217D-01 ,  2.368502055548D-01 ,  3.726715786250D-01 ,
c     * 5.092083225550D-01 ,  6.467110475491D-01 ,  7.854519690232D-01 ,
c     * 9.257184440339D-01 ,  1.067815135684D+00 ,  1.212067715646D+00 ,
c     * 1.358827381029D+00 ,  1.508476196330D+00 ,  1.661433435563D+00 ,
c     * 1.818163187438D+00 ,  1.979183577131D+00 ,  2.145078075547D+00 ,
c     * 2.316509528323D+00 ,  2.494237764181D+00 ,  2.679141968116D+00 ,
c     * 2.872249479966D+00 ,  3.074773384330D+00 ,  3.288162326323D+00 ,
c     * 3.514167642788D+00 ,  3.754935526960D+00 ,  4.013136236755D+00 ,
c     * 4.292149591579D+00 ,  4.596338646060D+00 ,  4.931466465190D+00 ,
c     * 5.305354988656D+00 ,  5.728974249695D+00 ,  6.218344086884D+00 ,
c     * 6.798089006814D+00 ,  7.508695084046D+00 ,  8.423188762282D+00 ,
c     * 9.692562072805D+00 ,  1.170765394670D+01 ,  1.607328906573D+01 /
C
C
C-----------------------------------------------------------------------
C
C
C
C
C                    Electron number density
      NSUBE = BRYDNS*YE
C
C                    Coefficants for chemical potential
C                    and thermodynamics quantities
      QSUBE = 1.0D0/( 3.0D0*(PI**2)*((HBAR*C)**3) )
C
      ACOEF = 0.5D0*NSUBE/QSUBE
C
      BCOEF = (ACOEF**2+((PI**6)*T**6)/27.0D0)**0.5D0
C
      DBDT = (PI**6)*(T**5)/(9.0D0*BCOEF)
C
      CCOEF = (ACOEF+BCOEF)**OVR3
C
C
C                    Electron chemical potential
      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF
C
C
C
C                            Initial guess at electron degeneracy
C                            parameter
      ETA_E = MUSUBE/T
C
      BETA = MASS_E/T
      DBTDT = -BETA/T
C                            Multiplicative factor in Fermi integrals
      TFAC = (T**3)/((PI**2)*(HBAR*C)**3)
      DTFDT = 3.0D0*TFAC/T
      PFAC = TFAC*T*OVR3
      EFAC = TFAC*T/BRYDNS
C
C                            Loop until the N-R iteration for the
C                            chemical potential converges
      DO 20 J=1,MAXITR,1
C
C                            Zero the electron number density integral
        NEN = 0.0D0
        DNENDE = 0.0D0
        DNENDB = 0.0D0
C
C                            Zero the electron energy density integral
        E_ENG = 0.0D0
        DEDE = 0.0D0
        DEDB = 0.0D0
C
C                            Zero the electron pressure integral
        E_PR = 0.0D0
        DPDE = 0.0D0
        DPDB = 0.0D0
C
C                            Zero the positron number density
        NEPLUS = 0.0D0
C
C       ----------------------------------------------------------------
C       |                    Do the Gaussian Quadrature                |
        DO 10 I=1,NGAUSS,1
C       |                    Electron occupation numbers               |
          FE = 1.0D0/(DEXP(X(I)+BETA-ETA_E)+1.0D0)
C       |                    Positron occupation numbers               |
          FP = 1.0D0/(DEXP(X(I)+BETA+ETA_E)+1.0D0)
C       |                    Derivatives w.r.t. beta                   |
          DFEDB = -FE*(1.0D0-FE)
          DFPDB = -FP*(1.0D0-FP)
C       |                    Derivatives w.r.t. Eta                    |
          DFEDE = FE*(1.0D0-FE)
          DFPDE = -FP*(1.0D0-FP)
C       |                                                              |
C       |                    Number density integral & derivatives     |
          XROOT = DSQRT(X(I)*(X(I)+2.0D0*BETA))
          XN1 = XROOT*(X(I)+BETA)
          DXN1DB = XROOT+X(I)*(X(I)+BETA)/XROOT
          NEN = NEN+WGHT(I)*XN1*(FE-FP)
          DNENDB = DNENDB+WGHT(I)*(DXN1DB*(FE-FP)+XN1*(DFEDB-DFPDB))
          DNENDE = DNENDE+WGHT(I)*XN1*(DFEDE-DFPDE)
C       |                    Sum for positron # integral               |
          NEPLUS = NEPLUS+WGHT(I)*XN1*FP
C       |                                                              |
C       |                                                              |
C       |                    Energy integral & derivatives             |
          XN2 = XROOT*(X(I)+BETA)**2
          DXN2DB = 2.0D0*XROOT*(X(I)+BETA)+
     1          X(I)*((X(I)+BETA)**2)/XROOT
          E_ENG = E_ENG+WGHT(I)*XN2*(FE+FP)
          DEDB = DEDB+WGHT(I)*(DXN2DB*(FP+FE)+XN2*(DFEDB+DFPDB))
          DEDE = DEDE+WGHT(I)*XN2*(DFEDE+DFPDE)
C       |                                                              |
C       |                                                              |
C       |                    Pressure integral & derivatives           |
          XN3 = XROOT**3
          DXN3DB = 3.0D0*XROOT*X(I)
          E_PR = E_PR+WGHT(I)*XN3*(FE+FP)
          DPDB = DPDB+WGHT(I)*(DXN3DB*(FE+FP)+XN3*(DFEDB+DFPDB))
          DPDE = DPDE+WGHT(I)*XN3*(DFEDE+DFPDE)
C       |                                                              |
C       |                                                              |
 10     CONTINUE
C       |                                                              |
C       ----------------------------------------------------------------
C                            Multiply by the temperature factor to get
C                            number density
        NE_CHK = TFAC*NEN
C                            Multiply by the temperature factor to get
C                            number density
        NEPLUS = TFAC*NEPLUS
C
C                            Calculate the new change in ETA
        D_ETA = -(NE_CHK-NSUBE)/(DNENDE*TFAC)
C
C                   If we've met the convergence criterion...
        IF(DABS(D_ETA).LT.(EPSIL*ETA_E)) THEN
C                          Then break out the chemical potential
C                          loop
          GOTO 30
        ELSE
C                          Otherwise update the chemical potential
          ETA_E = ETA_E+D_ETA
        ENDIF
C
 20   CONTINUE
C
C                                      If we reached this point the
C                                      N-R iteration didn't converge
      WRITE(*,*) ' EL_EOS: N-R iteration didnt converge!, 
     * T,YE,BRYDNS=',T,YE,BRYDNS
      STOP
C
C               If we reached this point then the N-R iteration has
C               converged.
 30   CONTINUE
C
C                     Is the result consistent with electron number
C                     density?
      IF(DABS((NE_CHK-NSUBE)/NSUBE).GT.1.0D-5) THEN
        WRITE(*,*) ' EL_EOS: Gaussian quadrature converged badly!'
        STOP
      ENDIF
C
C             Calculate thermodynamic quantities...
C
C                          Electron chemical potential
      MUSUBE = T*ETA_E
C                          Electron pressure
      EPRESS = PFAC*E_PR
C                          Electron internal energy per baryon
      EU = EFAC*E_ENG
C                          Electron free energy per baryon
      FSUBE = YE*MUSUBE-EPRESS/BRYDNS
C                          Electron entropy per baryon
      ES = (EU-FSUBE)/T
C
C                          Derivative of the electron eta w.r.t. T
      DETDT = -(DTFDT*NEN+TFAC*DNENDB*DBTDT)/(TFAC*DNENDE)
C                          Derivative of the electron eta w.r.t. NSUBE
      DETDNE = 1.0D0/(TFAC*DNENDE)
C
C                    Derivatives of chem. potential w.r.t. T,
C                    BRYDNS, YE
      DEMUDT = T*DETDT+MUSUBE/T
      DEMUDN = YE*T*DETDNE
      DEMUDY = BRYDNS*T*DETDNE
C
C
C                    Derivatives of pressure w.r.t. BRYDNS,YE,T
      DEPDN = YE*PFAC*DPDE*DETDNE
      DEPDY = BRYDNS*PFAC*DPDE*DETDNE
      DEPDT = (4.0D0*PFAC/T)*E_PR+PFAC*(DPDE*DETDT+DPDB*DBTDT)
C
C                    Derivatives of internal energy w.r.t.
C                    T,BRYDNS,YE
      DEUDT = (4.0D0*EFAC/T)*E_ENG+EFAC*(DEDB*DBTDT+DEDE*DETDT)
      DEUDN = (-1.0D0*EFAC/BRYDNS)*E_ENG+YE*EFAC*DEDE*DETDNE
      DEUDY = BRYDNS*EFAC*DEDE*DETDNE
C
C                    Derivatives of entropy w.r.t. T,BRYDNS,YE
      DESDT = -ES/T+(DEUDT-YE*DEMUDT+DEPDT/BRYDNS)/T
      DESDN = (DEUDN-YE*DEMUDN+DEPDN/BRYDNS-EPRESS/(BRYDNS**2))/T
      DESDY = (DEUDY-MUSUBE-YE*DEMUDY+DEPDY/BRYDNS)/T
C
C
C
C                    Photon pressure
      PPRESS = (PI**2)*(T**4)/(45.0*((HBAR*C)**3))
C                    Photon entropy per baryon
      PS = 4.0*PPRESS/(T*BRYDNS)
C
C                    Photon internal energy per baryon
      PU = 3.0*PPRESS/BRYDNS
C
C                    Photon free energy per baryon
      PF = PU-T*PS
C
C
C                    Derivatives of photon pressure
      DPPDN = 0.0
      DPPDT = BRYDNS*PS
      DPPDY = 0.0
C
C                    Derivatives of photon entropy
      DPSDN = -PS/BRYDNS
      DPSDT = 3.0*PS/T
      DPSDY = 0.0
C
C                    Derivatives of internal energy
      DPUDN = -0.75*T*PS/BRYDNS
      DPUDT = 3.0*PS
      DPUDY = 0.0
C23456789012345678901234567890123456789012345678901234567890123456789012
C
 999  RETURN
C
C
      END
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       F_2
C    TYPE:         DOUBLE PRECISION FUNCTION
C    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook
C    DATE:         12/16/91
C                  Version 2: 1/24/93
C
C    CALL LINE:    F_2(Y)      (2nd Fermi Integral)
C
C    INPUTS:       Y (DOUBLE PRECISION)   (Argument)
C
C    RETURN:       2nd Fermi Integral (DOUBLE PRECISION)
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION F_2(Y)
      IMPLICIT NONE
      DOUBLE PRECISION Y, YEXP
      IF(Y.GT.3.0D0) THEN
        WRITE(*,*) ' F_2(Y) FAILS FOR Y .GT. 3; Y =',Y
        STOP
      ENDIF
C
C                       Note: This approximation is based on the
C                       Bludman & Van Riper approximation (see
C                       Ap. J. Vol. 212 page 866-867 (1977))
C                       equation (3.6)
C
      IF(Y.LT.-1.0D0) THEN
        YEXP = EXP(Y)
        F_2 = 2.0*YEXP*(1.0-0.125*YEXP+0.037037*(YEXP**2))
      ELSE
        F_2 = 1.803D0+1.645D0*Y+0.6931D0*(Y**2)+0.1666667D0*(Y**3)+
     1        2.0833333D-2*(Y**4)-3.4722D-4*(Y**6)
      ENDIF
 999  RETURN
      END
 
 

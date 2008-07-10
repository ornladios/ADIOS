C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       EOS_M4C.INC
C    TYPE:         INCLUDE FILE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         6/27/90
C
C
C    CALL LINE:    INCLUDE 'EOS_M4C.INC/LIST'
C
C*************************************************************************
C
C                   INPUT & OUTPUT VARIBLES
C
C                            Electron fraction
      DOUBLE PRECISION YE
C                            Baryon density
      DOUBLE PRECISION BRYDNS
C                            Previous value of the nuclear proton
C                            fraction (must be supplied)
      DOUBLE PRECISION XPREV
C
C                            Previous value of the proton number
C                            density (must be supplied)
      DOUBLE PRECISION P_PREV
C
C                            Temperature, density, neutron chem. pot.,
C                            and proton chem. pot.
      DOUBLE PRECISION INPVAR(4)
C
C                            Flag to tell code what input is being
C                            supplied (i.e. what INPVAR is)
C                            1 --> Temperature
C                            2 --> Internal energy (NOT IMPLEMENTED)
C                            3 --> Entropy (NOT IMPLEMENTED)
      INTEGER IFLAG
C
C                            Flag returned by code to tell user what
C                            scheme was used
C                            1 --> "No Nuclei" scheme
C                            2 --> full scheme
C                            3 --> bulk scheme (above saturation dnsty)
      INTEGER EOSFLG
C
      INTEGER FORFLG, NGFLAG, NEWFLG, SSFLAG, FFLAG, ADFLAG
      INTEGER NF_FLG, RSFLAG, DBFLAG
cc      DATA NF_FLG/1/
cc      DATA RSFLAG/0/
C
C*************************************************************************
C                        Local variables
C*************************************************************************
C
C                            Flag used internally to this routine to
C                            record whether or not another scheme
C                            has already been tried
      INTEGER SWTFLG
C
C
C                            Accuracy in X
      DOUBLE PRECISION XACC
      PARAMETER (XACC=1.0D-8)
C
C                            Maximum number of X iterations
      INTEGER MAXIT, MAXIT2
      PARAMETER (MAXIT=200,MAXIT2=30)
C
C                            Update flag
      INTEGER UPFLAG
      PARAMETER(UPFLAG=3)
C
C                            External Interaction flag
C                            1 = on, 0 = off
      DOUBLE PRECISION EIFLAG
      PARAMETER(EIFLAG=1.0D+00)
C
C                            Alpha particle enable flag
C                            1= alphas, 0 = no alphas
      INTEGER ALFLAG
      PARAMETER(ALFLAG=1)
C
C                            Translational energy scaling factor
C                            1 = on, 0 = off
      DOUBLE PRECISION TRSCAL
      PARAMETER(TRSCAL=1.0D+00)
C
C                            Surface and Coulomb scaling factor
C                            1.0 = on, 0 = off
      DOUBLE PRECISION CSSCAL
      PARAMETER(CSSCAL=1.0D+00)
C
C
C
C                            Cutoff fraction of nuclei (If XH is below
C                            this the "no nuclei" scheme is used)
      DOUBLE PRECISION HEAVCT
      PARAMETER(HEAVCT=1.0D-6)
C
C                            Iteration loop variable
      INTEGER I, J, L, IT_NUM
C
C                            Combination of Coulomb and Surface
C                            coefficients
      DOUBLE PRECISION ALPHA, BETA
      PARAMETER (ALPHA=224.4D+00, BETA=1.3927D+00)
C
C                            Surface tension (MeV/fm**2)
      DOUBLE PRECISION SIG_0
      PARAMETER(SIG_0=46.5D+00)
C                            Binding energy of alpha particles (MeV)
      DOUBLE PRECISION BALPHA
      PARAMETER (BALPHA=28.3D+00)
      DOUBLE PRECISION V_ALFA
      PARAMETER (V_ALFA=24.0D+00)
C
C                            Nuclear symmetry energy (MeV)
      DOUBLE PRECISION SV
      PARAMETER (SV=29.3D+00)
C
C                            Nuclear compression modulus (MeV)
      DOUBLE PRECISION COMPRS, K
      PARAMETER (COMPRS=370.0D+00)
C
C                            Nuclear level density (per MeV)
      DOUBLE PRECISION LEVDNS, AV
      PARAMETER (LEVDNS=0.0666667D+00)
C
C                            Approximate size of nuclei
      DOUBLE PRECISION AZERO
      PARAMETER (AZERO=60.0D+00)
C
C
C                           Plancks constant & speed of light
      DOUBLE PRECISION HBAR, C
      PARAMETER (HBAR=6.58217317D-22,C=2.997924581D+23)
C
C                           Pi and square root of Pi
      DOUBLE PRECISION PI, ROOTPI
      PARAMETER(PI=3.1415927D+00,ROOTPI=1.7724539D+00)
C
C                           Square root of two
      DOUBLE PRECISION ROOT2
      PARAMETER (ROOT2=1.4142136D+00)
C
C                           Parameters containing powers of 1/3
C
      DOUBLE PRECISION OVR3, OVR23, OVR43, OVR53, OVR29
      DOUBLE PRECISION MOVR3, M2OVR3, M4OVR3, M5OVR3
      DOUBLE PRECISION OVR49
      PARAMETER(OVR49=0.44444444D+00)
      PARAMETER (OVR3=0.33333333D+00, OVR23=0.66666666D+00,
     *              OVR43=1.33333333D+00)
      PARAMETER (OVR53=1.66666667D+00, OVR29=0.2222222222D+00)
      PARAMETER (MOVR3=-0.33333333D+00, M2OVR3=-0.66666666D+00)
      PARAMETER (M4OVR3=-1.33333333D+00, M5OVR3=-1.666666667D+00)
C
C                          Ratio of baryon density to saturation density
      DOUBLE PRECISION Y
C
      DOUBLE PRECISION TFINAL
C
C                          Quantum concentration & Fermi integral
C                          coefficent
      DOUBLE PRECISION NQ, MQ, LQ, KQ, UQ
C
C                          Estimate of X (used in "no nuclei" scheme)
      DOUBLE PRECISION XGUESS
C
C                          Cutoff used to determine how to calc the
C                          number density of heavy nuclei (used in
C                          the "no nuclei" scheme)
      DOUBLE PRECISION MHCUT, MHOVT
      PARAMETER(MHCUT=30.0D+00)
C
C                          Chem. pot., mass, # density, & quantum
C                          concentration of heavy nuclei (used in
C                          the "no nuclei" scheme)
      DOUBLE PRECISION MUHEAV, NUCQ, MASSHV, NHEAVY
C
C
C
C
C                          Bulk nuclear pressures
      DOUBLE PRECISION PSUBI, POVRNI
C
C                          Neutron & proton rest masses
      DOUBLE PRECISION MASSN, MASSP
      PARAMETER(MASSN=939.5731D+00,MASSP=938.2796D+00)
C
C                          Saturation density, nuclear density, the
C                          difference of the two and their ratio
      DOUBLE PRECISION NSUBI, N_I, DNI, NRATIO, N_IOLD
      DOUBLE PRECISION NIOLD
C
C
C                          Surface energy coefficants
      DOUBLE PRECISION N_S
      DOUBLE PRECISION Q, R_0, ZETA_0
C
C
C
C                          Coulomb liquid correction functions
C                          & derivatives
      DOUBLE PRECISION W, DWDX, DWDT, DWDTDX, DWDXDX, DWDTDT
      DOUBLE PRECISION TZERO, DTZDX, DTZDXX
C
C
C                          Coulomb functions and derivatives
      DOUBLE PRECISION DU, DMU, DUP, DMUP, DUPP, DMUPP
      DOUBLE PRECISION DUX, DMUX, DUPX, DMUPX
      DOUBLE PRECISION DUT, DMUT, DUPT, DMUPT, DUPPT
      DOUBLE PRECISION DUXT, DMUXT
      DOUBLE PRECISION DUXX, DMUXX, DUTT, DMUTT
C
      DOUBLE PRECISION SCRDU, SCRDUP, SCRDPP, SCRD, SCRD2
      DOUBLE PRECISION SCRDUT, SCRDUX, SCRDXX, SCRDT, SCRDX, SCRDXT
      DOUBLE PRECISION SCRDTT
      DOUBLE PRECISION SCRD2T, SCRD2X, SCRDPT, SCRDPX
      DOUBLE PRECISION UTILDE, DUTIL, DUTILP, OMEG_R
C
C
C                          Finite size energy coefficent
      DOUBLE PRECISION ZETA, DZDT, DZDX, DZDNI
C
C                          Equilibrium equation finite size terms
      DOUBLE PRECISION BN, BP, BSUBP
C
C                          Baryon densities
      DOUBLE PRECISION NNOUT, ALFDNS, NUCDNS, NOUT
C
C                          Baryon fractions
      DOUBLE PRECISION XALFA, XNUT, XPROT, XH, XHCHK
      DOUBLE PRECISION XALFA2, XNUT2, XPROT2
C
C                          Fraction of space occupied by nuclei
      DOUBLE PRECISION U_NUC, U_N, COMPU, RY
C                          Fraction of space available to remaining
C                          baryons  (1-U_NUC)
      DOUBLE PRECISION EXALFA, EXCLU
C
C
C                          Surface tension variables
      DOUBLE PRECISION COMPX, SIGMA, AHBN, OVRX4, DSIGDX, DSGRDX, SIGSGP
      DOUBLE PRECISION SIGSG2
C
C
C                          Critical temperature variables
      DOUBLE PRECISION TSUBC, TSC1, DTCDX, DTCDXX, TSC_12
C
C                          Surface & translational temperature variables
      DOUBLE PRECISION DHDX, DHDXX, HX, DHDTDX, CHOVA
      DOUBLE PRECISION H, HPRIM, HPPRIM
      DOUBLE PRECISION CAPH, CAPHP
C
C                          Neutron-Proton mass difference
      DOUBLE PRECISION DELTAM
      PARAMETER(DELTAM=-1.2935D+00)
ccccccc      PARAMETER(DELTAM=0.0D+00)
C
C                          Chemical potentials
      DOUBLE PRECISION MUHAT, MUN, MUNOVT, MUALFA, MUPROT, MUN_I
      DOUBLE PRECISION MUP_I, MUN_O, MUP_O, MUN_Z, MUP_Z
      DOUBLE PRECISION ETA_NO, ETA_PO, ETA_NI, ETA_PI, ETAMAX
      DOUBLE PRECISION ETP_HI, ETP_LO, ETN_HI, ETN_LO
      DOUBLE PRECISION DELTMN, DELTMP, DELTAP
C
      DOUBLE PRECISION DMNDPI, DMPDPI, DPIDNI
C
      DOUBLE PRECISION OMEGA
C
      DOUBLE PRECISION DETA_P, DETA_N, NETA_N, NETA_P
      DOUBLE PRECISION CHNG_N, CHNG_P
C
C
C
C                          N-N interaction variables
      DOUBLE PRECISION VNOUT, VPOUT, VNI, VPI
      DOUBLE PRECISION ZNO, ZPO, ZNI, ZPI
C
C                          Nucleon kinetic energy densities
      DOUBLE PRECISION TAU_PO, TAU_NO, TAU_NI, TAU_PI
      DOUBLE PRECISION F32_PO, F32_NO, F32_NI, F32_PI
      DOUBLE PRECISION FRAT_N, FRAT_P
C
C                          Newton-Raphson equation variables
      DOUBLE PRECISION G, DGDX, DGDPRT, GOLD
C
C                          Atomic weight & number
      DOUBLE PRECISION A, Z
C
C                          Nuclear radius & volume
      DOUBLE PRECISION RSUBN, VSUBN
C
C                          Surface, Coulomb, translational, and
C                          bulk free energies (per baryon)
      DOUBLE PRECISION FSUBS, FSUBC, FSUBSC, FSUBI, FTRANS, MUSUBT
      DOUBLE PRECISION DMUTDT
C
      DOUBLE PRECISION FTR, DFTRDT, DFTRDU, DFTRDX, DFTRDN
      DOUBLE PRECISION F_SC, DFSCDT, DFSCDU, DFSCDX, DFSCDN
C
      DOUBLE PRECISION E_TR, DETRDT, DETRDU, DETRDX, DETRDN
      DOUBLE PRECISION E_SC, DESCDT, DESCDU, DESCDX, DESCDN
C
      DOUBLE PRECISION S_TR, DSTRDT, DSTRDU, DSTRDX, DSTRDN
      DOUBLE PRECISION S_SC, DSSCDT, DSSCDU, DSSCDX, DSSCDN
C
      DOUBLE PRECISION S_NUC, S_OUT, S_ALFA
      DOUBLE PRECISION E_NUC, E_OUT, E_ALFA
C
      DOUBLE PRECISION S_DENS, E_DENS
C
      DOUBLE PRECISION DETPDX, DETNDX, DETPDN, DETNDN
      DOUBLE PRECISION DTPIDX, DTNIDX, DTPDNI, DTNDNI
C
      DOUBLE PRECISION DSIDT, DSIDX, DSIDN, DEIDT, DEIDX, DEIDN
C
      DOUBLE PRECISION DSODT, DSODEP, DSODEN, DEODT, DEODEP, DEODEN
C
      DOUBLE PRECISION DSADT, DSADEP, DSADEN, DEADT, DEADEP, DEADEN
C
      DOUBLE PRECISION DNPDEP, DNPDEN, DNNDEP, DNNDEN
      DOUBLE PRECISION DTPDEP, DTPDEN, DTNDEP, DTNDEN
C
C
C
C
C
C                            BARYON THERMODYNAMIC VARIABLES
C
C                          Pressures
      DOUBLE PRECISION BPRESS, BPROUT, BPRNUC, BPRALF
C
C                          Entropies (per baryon)
      DOUBLE PRECISION BSOUT, BSNUC, BS, BSALFA
      DOUBLE PRECISION DSFSDU, DSFSDX, DSFSDN, DSFSDT
C
C
C
      DOUBLE PRECISION DBSODT, DBSNDT, DBSADT, SFS_DT
C
C                          Internal energies (per baryon)
      DOUBLE PRECISION BUOUT, BUNUC, BU, BUALFA
C
C                          Helmholtz free energies (per baryon)
      DOUBLE PRECISION BFOUT, BFALFA, BFNUC, BFTOT
C
C                          Adiabatic index
      DOUBLE PRECISION GAM_S
C
C                            TOTAL ENERGY, ENTROPY, AND PRESSURE
C
      DOUBLE PRECISION UTOT, STOT, PTOT, FTOT
      DOUBLE PRECISION FBARY, PBARY, MUBARY, MU_MAT
C
C                            X VARIABLES
C
C                          Internal proton fraction (in nuclei)
      DOUBLE PRECISION X
C                          Outside proton fraction (nucleon vapor)
      DOUBLE PRECISION XOUT
C
      DOUBLE PRECISION XNEW, XOLD, DX
C                          Limits on X
      DOUBLE PRECISION X_MIN
CC      DOUBLE PRECISION XL, XHI
CC      PARAMETER(XL=0.01D+00,XHI=0.52D+00)
C
C
C                            PROTON VARIABLES
      DOUBLE PRECISION NPOUT, PRTNEW, PRTOLD, DPRT
C
C
C                            TEMPERATURE VARIABLES
C
      DOUBLE PRECISION T
C
C                            Diagonal line in the temp-density plane
      DOUBLE PRECISION TCHECK
C
C                            Temporary storage variables
      DOUBLE PRECISION TMP1, TMP1P, TMP1PP
      DOUBLE PRECISION TMP2, TMP2P, TMP2PP, TMP2PT, TMP2PX
      DOUBLE PRECISION TMP2T, TMP2X, TMP2XX, TMP2XT, TMP2TT
      DOUBLE PRECISION TMP3, TMP3P, TMP3PP, TMP3PT, TMP3PX
      DOUBLE PRECISION TMP3T, TMP3X, TMP3XX, TMP3XT, TMP3TT
      DOUBLE PRECISION TMP4, TMP5
C
C
C
      DOUBLE PRECISION NSUBIN, NETAP, NETAN
C
      DOUBLE PRECISION NSIACC, PRTACC, NUTACC
      PARAMETER(NSIACC=1.0D-8,PRTACC=1.0D-9,NUTACC=1.0D-9)
cccc      PARAMETER(NSIACC=1.0D-9,PRTACC=1.0D-10,NUTACC=1.0D-10)
C
      DOUBLE PRECISION DETERM
C
      DOUBLE PRECISION DNSUBI, DETAP, DETAN
C
      DOUBLE PRECISION A1, A2, A3
C
      DOUBLE PRECISION A1H, A1L, DA1DN, DA1ETP, DA1ETN
      DOUBLE PRECISION A2H, A2L, DA2DN, DA2ETP, DA2ETN
      DOUBLE PRECISION A3H, A3L, DA3DN, DA3ETP, DA3ETN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION GPI,GNI,GPO,GNO
      DOUBLE PRECISION DVPIDP, DVPIDN, DVNIDP, DVNIDN
      DOUBLE PRECISION DVPODP, DVPODN, DVNODP, DVNODN
      DOUBLE PRECISION MSSCON
      DOUBLE PRECISION DTPIDP, DTPIDN, DTNIDP, DTNIDN
      DOUBLE PRECISION DTPODP, DTPODN, DTNODP, DTNODN
      DOUBLE PRECISION DMPIDP, DMPIDN, DMNIDP, DMNIDN
      DOUBLE PRECISION DMPODP, DMPODN, DMNODP, DMNODN
      DOUBLE PRECISION DPIDP, DPIDN, DPODP, DPODN, DPADP, DPADN
      DOUBLE PRECISION N1, N2, DUDPO, DUDNO, DUDNI
      DOUBLE PRECISION DXDPO, DXDNO, DXDNI
      DOUBLE PRECISION DB1DNI, DB1DX, DB1DU
      DOUBLE PRECISION DB2DNI, DB2DX, DB2DU
      DOUBLE PRECISION DB3DNI, DB3DX, DB3DU
      DOUBLE PRECISION DA1ID1, DA1ID2, DA1ID3, DA1OD1, DA1OD2, DA1OD3
      DOUBLE PRECISION DA2ID1, DA2ID2, DA2ID3, DA2OD1, DA2OD2, DA2OD3
      DOUBLE PRECISION DA3ID1, DA3ID2, DA3ID3, DA3OD1, DA3OD2, DA3OD3
      DOUBLE PRECISION DA1D1, DA1D2, DA1D3, DB1D1, DB1D2, DB1D3
      DOUBLE PRECISION DA2D1, DA2D2, DA2D3, DB2D1, DB2D2, DB2D3
      DOUBLE PRECISION DA3D1, DA3D2, DA3D3, DB3D1, DB3D2, DB3D3
      DOUBLE PRECISION DNDETN, DPDETP
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
C/////////////////////////////////////////////////////////////////////
C                   Derivative variables
C/////////////////////////////////////////////////////////////////////
C
C
C                   Temperature derivatives
C
      DOUBLE PRECISION DMPODT, DMNODT, DMPIDT, DMNIDT
      DOUBLE PRECISION DNADPO, DNADNO
C
      DOUBLE PRECISION DTPIDT, DTNIDT, DTPODT, DTNODT
      DOUBLE PRECISION DPODT, DPIDT
C
      DOUBLE PRECISION DNPODT, DNNODT
C
      DOUBLE PRECISION DV_DT
C
      DOUBLE PRECISION DNADT, DPADT, DMUADT
C
C
      DOUBLE PRECISION DU_DT, DX_DT, DNI_DT, DEP_DT, DEN_DT
C
      DOUBLE PRECISION DU_DN, DX_DN, DNI_DN, DEP_DN, DEN_DN
C
      DOUBLE PRECISION DU_DY, DX_DY, DNI_DY, DEP_DY, DEN_DY
C
      DOUBLE PRECISION DNA_DT, DNA_DY, DNA_DN
C
      DOUBLE PRECISION DB1DT, DB2DT, DB3DT
C
C
C
C                   Exterior particle density derivatives
C
C
      DOUBLE PRECISION DV_DPO, DV_DNO
C
C
C
C
C                   ETA derivatives
C
C
      DOUBLE PRECISION DMPDEP, DMPDEN, DMNDEP, DMNDEN
C
      DOUBLE PRECISION DV_DEP, DV_DEN
C
      DOUBLE PRECISION DNADEP, DNADEN, DPADEP, DPADEN
C
      DOUBLE PRECISION DPODEP, DPODEN
C
      DOUBLE PRECISION DNPODN, DNPIDN, DNNIDN
C
      DOUBLE PRECISION DMADEP, DMADEN
C
C
C
      DOUBLE PRECISION DFDOM(5,5), DFDOMI(5,5)
      DOUBLE PRECISION DFDL_1(5), DFDL_2(5), DFDL_3(5)
C
      DOUBLE PRECISION DG1DO1, DG1DO2, DG2DO1, DG2DO2, DET_GT
      DOUBLE PRECISION DG1DL1, DG1DL2, DG1DL3
      DOUBLE PRECISION DG2DL1, DG2DL2, DG2DL3
C
      DOUBLE PRECISION A_1X, A_1U, A_1N
      DOUBLE PRECISION A_2X, A_2U, A_2N
      DOUBLE PRECISION A_3X, A_3U, A_3N
      DOUBLE PRECISION DET_AL
C
      DOUBLE PRECISION AT_11, AT_12, AT_13
      DOUBLE PRECISION AT_21, AT_22, AT_23
      DOUBLE PRECISION AT_31, AT_32, AT_33
C
      DOUBLE PRECISION B_1P, B_2P, B_3P
      DOUBLE PRECISION B_1N, B_2N, B_3N
C
C
      DOUBLE PRECISION DNIDPO, DNIDNO
C
C
      DOUBLE PRECISION C_P, C_N, CP_P, CP_N
      DOUBLE PRECISION BB_P, BB_N, BP_P, BP_N
C
      DOUBLE PRECISION DMP_DN, DMN_DN
C
      DOUBLE PRECISION DMP_DY, DMN_DY
C
      DOUBLE PRECISION DFDNDN, DFDYDY, DFDNDY, DFDNDT, DFDYDT
C
      DOUBLE PRECISION DFDTDT
C
      DOUBLE PRECISION DBPDT, DBPDN, DBPDY
C
      DOUBLE PRECISION DBSDT, DBSDN, DBSDY
C
      DOUBLE PRECISION DBUDT, DBUDN, DBUDY
C
      DOUBLE PRECISION DBMUDT, DBMUDN, DBMUDY
C
      DOUBLE PRECISION DMUDT, DMUDN, DMUDY, DSDT, DSDN, DSDY
      DOUBLE PRECISION DPDT, DPDN, DPDY, DUDT, DUDN, DUDY
C
C
      DOUBLE PRECISION NPI, NNI, DNPIDT, DNNIDT
C
      DOUBLE PRECISION DXHDT, DXODT, DXADT
C
      DOUBLE PRECISION VP_IN, VN_IN
C
C
      DOUBLE PRECISION DUODT, DUIDT, DUADT
C
      DOUBLE PRECISION DBFDT, DBFDY, DBFDN
C
C
C
C
C
C
C/////////////////////////////////////////////////////////////////////
C
C
C
C
C                   This common block contains all of the thermodynamic
C                   and compositional variables that the EOS code
C                   calculates & it should be included anywhere
C                   these variables need to be referenced
C
C
      COMMON /EOSVAR/ MUN,MUN_I,MUPROT,MUALFA,MUHAT,
     1    BU, BS, BPRESS, UTOT, STOT, PTOT,
     2    XPROT, XNUT, XALFA, XH, X, A,
     3    BUOUT, BUNUC, BPROUT, BPRNUC, BPRALF, BSOUT,BSNUC,BFOUT,
     4    BFNUC, BFALFA, FTOT, BFTOT,U_NUC, NSUBI, DNI, RSUBN,
     5    VSUBN, FBARY, PBARY, MUBARY, MU_MAT,
     6    XALFA2, XNUT2, XPROT2, IT_NUM
C
C
C                    This common contains the variables that are used
C                    in solving the chemical and pressure equilibrium
C                    equations
C
      COMMON /EQLEQN/ NNOUT, NPOUT, NOUT, VNOUT, VPOUT, F32_NO, F32_PO,
     1 MUN_O, MUP_O, ALFDNS, COMPX, SIGMA, SIGSGP, DU, DMU, DUP, DMUP,
     2 SCRDU, SCRDUP, SCRD, ZETA, FSUBS, FSUBC, MUSUBT, FTRANS, NRATIO,
     3 VNI, VPI, ETA_NI, ETA_PI, MUP_I, F32_NI, F32_PI, BN, BP, BSUBP,
     4 H, HPRIM, PSUBI, TSUBC, DHDX, DHDXX, SCRD2, SCRDPP, TFINAL
C
C
C                   This common block contains the derivatives of
C                   the thermodynamic potentials (both baryon and
C                   the totals)
C
      COMMON /DERIVS/ DBMUDT, DBMUDN, DBMUDY, DMUDT, DMUDN, DMUDY,
     1                DBPDT, DBPDN, DBPDY, DPDT, DPDN, DPDY,
     2                DBSDT, DBSDN, DBSDY, DSDT, DSDN, DSDY,
     3                DBUDT, DBUDN, DBUDY, DUDT, DUDN, DUDY, GAM_S,
     4                DU_DT, DX_DT, DNI_DT, DEP_DT, DEN_DT,
     5                DU_DY, DX_DY, DNI_DY, DEP_DY, DEN_DY,
     6                DU_DN, DX_DN, DNI_DN, DEP_DN, DEN_DN
C
C
C                   This common block contains control flags
      COMMON /FLAGS/ ADFLAG, NF_FLG, RSFLAG, DBFLAG
C
C
C                   This common block contains the quality of solution
C                   variables that determine how "zeroed" the equilbrium
C                   equations are
      COMMON /QCHECK/ A1,A2,A3
C
      double precision ncomp
      COMMON /DTEST/ DFDOM, DFDL_1,DFDL_2,DFDL_3,DNA_DT,DNA_DN,DNA_DY,
     1 ncomp
C

!-----------------------------------------------------------------------
!    MODULE:       eos_m4c_module
!    Author:       S. W. Bruenn
!    Date:         10/26/02
!-----------------------------------------------------------------------
MODULE eos_m4c_module

USE kind_module
SAVE

!-----------------------------------------------------------------------
!
!               \\\\\ INPUT & OUTPUT VARIBLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Electron fraction
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: YE

!-----------------------------------------------------------------------
!  Baryon density
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: BRYDNS

!-----------------------------------------------------------------------
!  Previous value of the nuclear proton fraction (must be supplied)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: XPREV

!-----------------------------------------------------------------------
!  Previous value of the proton number density (must be supplied)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: P_PREV

!-----------------------------------------------------------------------
!  Temperature, density, neutron chem. pot., and proton chem. pot
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(4)                 :: INPVAR

!-----------------------------------------------------------------------
!  Flag to tell code what input is being supplied (i.e. what INPVAR is)
!     1 --> Temperature
!     2 --> Internal energy (NOT IMPLEMENTED)
!     3 --> Entropy (NOT IMPLEMENTED)
!-----------------------------------------------------------------------

INTEGER                                         :: IFLAG

!-----------------------------------------------------------------------
!  Flag returned by code to tell user what scheme was used
!     1 --> "No Nuclei" scheme
!     2 --> full scheme
!     3 --> bulk scheme (above saturation dnsty)
!-----------------------------------------------------------------------

INTEGER                                         :: EOSFLG

INTEGER                                         :: FORFLG, NGFLAG, NEWFLG, SSFLAG, FFLAG, ADFLAG
INTEGER                                         :: NF_FLG, RSFLAG, DBFLAG

!-----------------------------------------------------------------------
!
!                    \\\\\ LOCAL VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Flag used internally to this routine to record whether or not another
!   scheme has already been tried
!-----------------------------------------------------------------------

INTEGER                                         :: SWTFLG

!-----------------------------------------------------------------------
!  XACC - Accuracy in X
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: XACC = 1.0D-8

!-----------------------------------------------------------------------
!  Maximum number of X iterations
!-----------------------------------------------------------------------

INTEGER, PARAMETER                              :: MAXIT = 200

!-----------------------------------------------------------------------
!  Update flag
!-----------------------------------------------------------------------

INTEGER, PARAMETER                              :: UPFLAG = 3

!-----------------------------------------------------------------------
!  External Interaction flag
!     1 = on, 0 = off
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: EIFLAG = 1.0D+00

!-----------------------------------------------------------------------
!  Alpha particle enable flag
!     1= alphas, 0 = no alphas
!-----------------------------------------------------------------------

INTEGER, PARAMETER                              :: ALFLAG = 1

!-----------------------------------------------------------------------
!  Translational energy scaling factor
!     1 = on, 0 = off
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: TRSCAL = 1.0D+00

!-----------------------------------------------------------------------
!  Surface and Coulomb scaling factor
!     1.0 = on, 0 = off
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: CSSCAL = 1.0D+00

!-----------------------------------------------------------------------
!  Cutoff fraction of nuclei (If XH is below this the "no nuclei" scheme
!   is used)
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: HEAVCT = 1.0D-6

!-----------------------------------------------------------------------
!  Iteration loop variables
!-----------------------------------------------------------------------

INTEGER                                         :: I, J, L, IT_NUM

!-----------------------------------------------------------------------
!  Combination of Coulomb and Surface coefficients
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: ALPHA = 224.4D+00
REAL(KIND=double), PARAMETER                    :: BETA = 1.3927D+00

!-----------------------------------------------------------------------
!  Surface tension (MeV/fm**2)
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: SIG_0 = 46.5D+00

!-----------------------------------------------------------------------
!  Binding energy of alpha particles (MeV)
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: BALPHA = 28.3D+00
REAL(KIND=double), PARAMETER                    :: V_ALFA = 24.0D+00

!-----------------------------------------------------------------------
!  Nuclear symmetry energy (MeV)
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: SV = 29.3D+00

!-----------------------------------------------------------------------
!  Nuclear compression modulus (MeV)
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: COMPRS = 370.0D+00
REAL(KIND=double)                               :: K

!-----------------------------------------------------------------------
!  Nuclear level density (per MeV)
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: LEVDNS = 0.0666667D+00
REAL(KIND=double)                               :: AV

!-----------------------------------------------------------------------
!  Approximate size of nuclei
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: AZERO = 60.0D+00

!-----------------------------------------------------------------------
!  Plancks constant & speed of light
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: HBAR = 6.58217317D-22
REAL(KIND=double), PARAMETER                    :: C = 2.997924581D+23

!-----------------------------------------------------------------------
!  Pi and square root of Pi
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: PI = 3.1415927D+00
REAL(KIND=double), PARAMETER                    :: ROOTPI = 1.7724539D+00

!-----------------------------------------------------------------------
!  Square root of two
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: ROOT2 = 1.4142136D+00

!-----------------------------------------------------------------------
!  Parameters containing powers of 1/3
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: OVR3 = 0.33333333D+00
REAL(KIND=double), PARAMETER                    :: OVR23 = 0.66666666D+00
REAL(KIND=double), PARAMETER                    :: OVR43 = 1.33333333D+00
REAL(KIND=double), PARAMETER                    :: OVR53 = 1.66666667D+00
REAL(KIND=double), PARAMETER                    :: OVR29 = 0.2222222222D+00
REAL(KIND=double), PARAMETER                    :: OVR49 = 0.44444444D+00
REAL(KIND=double), PARAMETER                    :: MOVR3 = -0.33333333D+00
REAL(KIND=double), PARAMETER                    :: M2OVR3 = -0.66666666D+00
REAL(KIND=double), PARAMETER                    :: M4OVR3 = -1.33333333D+00
REAL(KIND=double), PARAMETER                    :: M5OVR3 = -1.666666667D+00

!-----------------------------------------------------------------------
!  Ratio of baryon density to saturation density
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: Y
REAL(KIND=double)                               :: TFINAL

!-----------------------------------------------------------------------
!  Quantum concentration & Fermi integral coefficent
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: NQ, MQ, LQ, KQ, UQ

!-----------------------------------------------------------------------
!  Estimate of X (used in "no nuclei" scheme)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: XGUESS

!-----------------------------------------------------------------------
!  Cutoff used to determine how to calc the number density of heavy
!   nuclei (used in the "no nuclei" scheme)
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: MHCUT = 30.0D+00
REAL(KIND=double)                               :: MHOVT

!-----------------------------------------------------------------------
!  Chem. pot., mass, # density, & quantum concentration of heavy nuclei
!   (used in the "no nuclei" scheme)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: MUHEAV, NUCQ, MASSHV, NHEAVY

!-----------------------------------------------------------------------
!  Bulk nuclear pressures
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: PSUBI, POVRNI

!-----------------------------------------------------------------------
!  Neutron & proton rest masses
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: MASSN = 939.5731D+00
REAL(KIND=double), PARAMETER                    :: MASSP = 938.2796D+00

!-----------------------------------------------------------------------
!  Saturation density, nuclear density, the difference of the two and
!   their ratio
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: NSUBI, N_I, DNI, NRATIO, N_IOLD, NIOLD

!-----------------------------------------------------------------------
!  Surface energy coefficants
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: N_S, Q, R_0, ZETA_0

!-----------------------------------------------------------------------
!  Coulomb liquid correction functions & derivatives
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: W, DWDX, DWDT, DWDTDX, DWDXDX, DWDTDT
REAL(KIND=double)                               :: TZERO, DTZDX, DTZDXX

!-----------------------------------------------------------------------
!  Coulomb functions and derivatives
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: DU, DMU, DUP, DMUP, DUPP, DMUPP
REAL(KIND=double)                               :: DUX, DMUX, DUPX, DMUPX
REAL(KIND=double)                               :: DUT, DMUT, DUPT, DMUPT, DUPPT
REAL(KIND=double)                               :: DUXT, DMUXT
REAL(KIND=double)                               :: DUXX, DMUXX, DUTT, DMUTT

REAL(KIND=double)                               :: SCRDU, SCRDUP, SCRDPP, SCRD, SCRD2
REAL(KIND=double)                               :: SCRDUT, SCRDUX, SCRDXX, SCRDT, SCRDX, SCRDXT
REAL(KIND=double)                               :: SCRDTT
REAL(KIND=double)                               :: SCRD2T, SCRD2X, SCRDPT, SCRDPX
REAL(KIND=double)                               :: UTILDE, DUTIL, DUTILP, OMEG_R

!-----------------------------------------------------------------------
!  Finite size energy coefficent
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: ZETA, DZDT, DZDX, DZDNI

!-----------------------------------------------------------------------
!  Equilibrium equation finite size terms
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: BN, BP, BSUBP

!-----------------------------------------------------------------------
!  Baryon densities
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: NNOUT, ALFDNS, NUCDNS, NOUT

!-----------------------------------------------------------------------
!  Baryon fractions
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: XALFA, XNUT, XPROT, XH, XHCHK
REAL(KIND=double)                               :: XALFA2, XNUT2, XPROT2

!-----------------------------------------------------------------------
!  Fraction of space occupied by nuclei
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: U_NUC, U_N, COMPU, RY

!-----------------------------------------------------------------------
!  Fraction of space available to remaining baryons  (1-U_NUC)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: EXALFA, EXCLU

!-----------------------------------------------------------------------
!  Surface tension variables
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: COMPX, SIGMA, AHBN, OVRX4, DSIGDX, DSGRDX, SIGSGP
REAL(KIND=double)                               :: SIGSG2

!-----------------------------------------------------------------------
!  Critical temperature variables
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: TSUBC, TSC1, DTCDX, DTCDXX, TSC_12

!-----------------------------------------------------------------------
!  Surface & translational temperature variables
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: DHDX, DHDXX, HX, DHDTDX, CHOVA
REAL(KIND=double)                               :: H, HPRIM, HPPRIM
REAL(KIND=double)                               :: CAPH, CAPHP

!-----------------------------------------------------------------------
!  Neutron-Proton mass difference
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER                    :: DELTAM = -1.2935D+00

!-----------------------------------------------------------------------
!  Chemical potentials
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: MUHAT, MUN, MUNOVT, MUALFA, MUPROT, MUN_I
REAL(KIND=double)                               :: MUP_I, MUN_O, MUP_O, MUN_Z, MUP_Z
REAL(KIND=double)                               :: ETA_NO, ETA_PO, ETA_NI, ETA_PI, ETAMAX
REAL(KIND=double)                               :: ETP_HI, ETP_LO, ETN_HI, ETN_LO
REAL(KIND=double)                               :: DELTMN, DELTMP, DELTAP

REAL(KIND=double)                               :: DMNDPI, DMPDPI, DPIDNI

REAL(KIND=double)                               :: OMEGA

REAL(KIND=double)                               :: DETA_P, DETA_N, NETA_N, NETA_P
REAL(KIND=double)                               :: CHNG_N, CHNG_P

!-----------------------------------------------------------------------
!  Chemical potentials
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: VNOUT, VPOUT, VNI, VPI
REAL(KIND=double)                               :: ZNO, ZPO, ZNI, ZPI

!-----------------------------------------------------------------------
!  Nucleon kinetic energy densities
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: TAU_PO, TAU_NO, TAU_NI, TAU_PI
REAL(KIND=double)                               :: F32_PO, F32_NO, F32_NI, F32_PI
REAL(KIND=double)                               :: FRAT_N, FRAT_P

!-----------------------------------------------------------------------
!  Newton-Raphson equation variables
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: G, DGDX, DGDPRT, GOLD

!-----------------------------------------------------------------------
!  Atomic weight & number
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: A, Z

!-----------------------------------------------------------------------
!  Nuclear radius & volume
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: RSUBN, VSUBN

!-----------------------------------------------------------------------
!  Surface, Coulomb, translational, and. bulk free energies (per baryon)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: FSUBS, FSUBC, FSUBSC, FSUBI, FTRANS, MUSUBT
REAL(KIND=double)                               :: DMUTDT

REAL(KIND=double)                               :: FTR, DFTRDT, DFTRDU, DFTRDX, DFTRDN
REAL(KIND=double)                               :: F_SC, DFSCDT, DFSCDU, DFSCDX, DFSCDN

REAL(KIND=double)                               :: E_TR, DETRDT, DETRDU, DETRDX, DETRDN
REAL(KIND=double)                               :: E_SC, DESCDT, DESCDU, DESCDX, DESCDN

REAL(KIND=double)                               :: S_TR, DSTRDT, DSTRDU, DSTRDX, DSTRDN
REAL(KIND=double)                               :: S_SC, DSSCDT, DSSCDU, DSSCDX, DSSCDN

REAL(KIND=double)                               :: S_NUC, S_OUT, S_ALFA
REAL(KIND=double)                               :: E_NUC, E_OUT, E_ALFA

REAL(KIND=double)                               :: S_DENS, E_DENS

REAL(KIND=double)                               :: DETPDX, DETNDX, DETPDN, DETNDN
REAL(KIND=double)                               :: DTPIDX, DTNIDX, DTPDNI, DTNDNI

REAL(KIND=double)                               :: DSIDT, DSIDX, DSIDN, DEIDT, DEIDX, DEIDN

REAL(KIND=double)                               :: DSODT, DSODEP, DSODEN, DEODT, DEODEP, DEODEN

REAL(KIND=double)                               :: DSADT, DSADEP, DSADEN, DEADT, DEADEP, DEADEN

REAL(KIND=double)                               :: DNPDEP, DNPDEN, DNNDEP, DNNDEN
REAL(KIND=double)                               :: DTPDEP, DTPDEN, DTNDEP, DTNDEN

!-----------------------------------------------------------------------
!
!               \\\\\ BARYON THERMODYNAMIC VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Pressures
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: BPRESS, BPROUT, BPRNUC, BPRALF

!-----------------------------------------------------------------------
!  Entropies (per baryon)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: BSOUT, BSNUC, BS, BSALFA
REAL(KIND=double)                               :: DSFSDU, DSFSDX, DSFSDN, DSFSDT

REAL(KIND=double)                               :: DBSODT, DBSNDT, DBSADT, SFS_DT

!-----------------------------------------------------------------------
!  Internal energies (per baryon)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: BUOUT, BUNUC, BU, BUALFA

!-----------------------------------------------------------------------
!  Helmholtz free energies (per baryon)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: BFOUT, BFALFA, BFNUC, BFTOT

!-----------------------------------------------------------------------
!  Adiabatic index
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: GAM_S

!-----------------------------------------------------------------------
!
!            \\\\\ TOTAL ENERGY, ENTROPY, AND PRESSURE /////
!
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: UTOT, STOT, PTOT, FTOT
REAL(KIND=double)                               :: FBARY, PBARY, MUBARY, MU_MAT

!-----------------------------------------------------------------------
!
!                       \\\\\ X VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Internal proton fraction (in nuclei)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: X

!-----------------------------------------------------------------------
!  Outside proton fraction (nucleon vapor)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: XOUT

REAL(KIND=double)                               :: XNEW, XOLD, DX

!-----------------------------------------------------------------------
!  Limits on X
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: X_MIN

!-----------------------------------------------------------------------
!
!                     \\\\\ PROTON VARIABLES /////
!
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: NPOUT, PRTNEW, PRTOLD, DPRT

!-----------------------------------------------------------------------
!
!                   \\\\\ TEMPERATURE VARIABLES /////
!
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: T

!-----------------------------------------------------------------------
!  Diagonal line in the temp-density plane
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: TCHECK

!-----------------------------------------------------------------------
!  Temporary storage variables
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: TMP1, TMP1P, TMP1PP
REAL(KIND=double)                               :: TMP2, TMP2P, TMP2PP, TMP2PT, TMP2PX
REAL(KIND=double)                               :: TMP2T, TMP2X, TMP2XX, TMP2XT, TMP2TT
REAL(KIND=double)                               :: TMP3, TMP3P, TMP3PP, TMP3PT, TMP3PX
REAL(KIND=double)                               :: TMP3T, TMP3X, TMP3XX, TMP3XT, TMP3TT
REAL(KIND=double)                               :: TMP4, TMP5

REAL(KIND=double)                               :: NSUBIN, NETAP, NETAN

REAL(KIND=double), PARAMETER                    :: NSIACC = 1.0D-8
REAL(KIND=double), PARAMETER                    :: PRTACC = 1.0D-9
REAL(KIND=double), PARAMETER                    :: NUTACC = 1.0D-9

REAL(KIND=double)                               :: DETERM

REAL(KIND=double)                               :: DNSUBI, DETAP, DETAN

REAL(KIND=double)                               :: A1, A2, A3

REAL(KIND=double)                               :: A1H, A1L, DA1DN, DA1ETP, DA1ETN
REAL(KIND=double)                               :: A2H, A2L, DA2DN, DA2ETP, DA2ETN
REAL(KIND=double)                               :: A3H, A3L, DA3DN, DA3ETP, DA3ETN
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
REAL(KIND=double)                               :: GPI,GNI,GPO,GNO
REAL(KIND=double)                               :: DVPIDP, DVPIDN, DVNIDP, DVNIDN
REAL(KIND=double)                               :: DVPODP, DVPODN, DVNODP, DVNODN
REAL(KIND=double)                               :: MSSCON
REAL(KIND=double)                               :: DTPIDP, DTPIDN, DTNIDP, DTNIDN
REAL(KIND=double)                               :: DTPODP, DTPODN, DTNODP, DTNODN
REAL(KIND=double)                               :: DMPIDP, DMPIDN, DMNIDP, DMNIDN
REAL(KIND=double)                               :: DMPODP, DMPODN, DMNODP, DMNODN
REAL(KIND=double)                               :: DPIDP, DPIDN, DPODP, DPODN, DPADP, DPADN
REAL(KIND=double)                               :: N1, N2, DUDPO, DUDNO, DUDNI
REAL(KIND=double)                               :: DXDPO, DXDNO, DXDNI
REAL(KIND=double)                               :: DB1DNI, DB1DX, DB1DU
REAL(KIND=double)                               :: DB2DNI, DB2DX, DB2DU
REAL(KIND=double)                               :: DB3DNI, DB3DX, DB3DU
REAL(KIND=double)                               :: DA1ID1, DA1ID2, DA1ID3, DA1OD1, DA1OD2, DA1OD3
REAL(KIND=double)                               :: DA2ID1, DA2ID2, DA2ID3, DA2OD1, DA2OD2, DA2OD3
REAL(KIND=double)                               :: DA3ID1, DA3ID2, DA3ID3, DA3OD1, DA3OD2, DA3OD3
REAL(KIND=double)                               :: DA1D1, DA1D2, DA1D3, DB1D1, DB1D2, DB1D3
REAL(KIND=double)                               :: DA2D1, DA2D2, DA2D3, DB2D1, DB2D2, DB2D3
REAL(KIND=double)                               :: DA3D1, DA3D2, DA3D3, DB3D1, DB3D2, DB3D3
REAL(KIND=double)                               :: DNDETN, DPDETP
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!
!/////////////////////////////////////////////////////////////////////
!                   Derivative variables
!/////////////////////////////////////////////////////////////////////
!
!

!-----------------------------------------------------------------------
!  Temperature derivatives
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: DMPODT, DMNODT, DMPIDT, DMNIDT
REAL(KIND=double)                               :: DNADPO, DNADNO

REAL(KIND=double)                               :: DTPIDT, DTNIDT, DTPODT, DTNODT
REAL(KIND=double)                               :: DPODT, DPIDT

REAL(KIND=double)                               :: DNPODT, DNNODT

REAL(KIND=double)                               :: DV_DT

REAL(KIND=double)                               :: DNADT, DPADT, DMUADT

REAL(KIND=double)                               :: DU_DT, DX_DT, DNI_DT, DEP_DT, DEN_DT

REAL(KIND=double)                               :: DU_DN, DX_DN, DNI_DN, DEP_DN, DEN_DN

REAL(KIND=double)                               :: DU_DY, DX_DY, DNI_DY, DEP_DY, DEN_DY

REAL(KIND=double)                               :: DNA_DT, DNA_DY, DNA_DN

REAL(KIND=double)                               :: DB1DT, DB2DT, DB3DT

!-----------------------------------------------------------------------
!  Exterior particle density derivatives
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: DV_DPO, DV_DNO


!-----------------------------------------------------------------------
!  ETA derivatives
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: DMPDEP, DMPDEN, DMNDEP, DMNDEN

REAL(KIND=double)                               :: DV_DEP, DV_DEN

REAL(KIND=double)                               :: DNADEP, DNADEN, DPADEP, DPADEN

REAL(KIND=double)                               :: DPODEP, DPODEN

REAL(KIND=double)                               :: DNPODN, DNPIDN, DNNIDN

REAL(KIND=double)                               :: DMADEP, DMADEN

REAL(KIND=double), DIMENSION(5,5)               :: DFDOM, DFDOMI
REAL(KIND=double), DIMENSION(5)                 :: DFDL_1, DFDL_2, DFDL_3

REAL(KIND=double)                               :: DG1DO1, DG1DO2, DG2DO1, DG2DO2, DET_GT
REAL(KIND=double)                               :: DG1DL1, DG1DL2, DG1DL3
REAL(KIND=double)                               :: DG2DL1, DG2DL2, DG2DL3

REAL(KIND=double)                               :: A_1X, A_1U, A_1N
REAL(KIND=double)                               :: A_2X, A_2U, A_2N
REAL(KIND=double)                               :: A_3X, A_3U, A_3N
REAL(KIND=double)                               :: DET_AL

REAL(KIND=double)                               :: AT_11, AT_12, AT_13
REAL(KIND=double)                               :: AT_21, AT_22, AT_23
REAL(KIND=double)                               :: AT_31, AT_32, AT_33

REAL(KIND=double)                               :: B_1P, B_2P, B_3P
REAL(KIND=double)                               :: B_1N, B_2N, B_3N

REAL(KIND=double)                               :: DNIDPO, DNIDNO

REAL(KIND=double)                               :: C_P, C_N, CP_P, CP_N
REAL(KIND=double)                               :: BB_P, BB_N, BP_P, BP_N

REAL(KIND=double)                               :: DMP_DN, DMN_DN

REAL(KIND=double)                               :: DMP_DY, DMN_DY

REAL(KIND=double)                               :: DFDNDN, DFDYDY, DFDNDY, DFDNDT, DFDYDT

REAL(KIND=double)                               :: DFDTDT

REAL(KIND=double)                               :: DBPDT, DBPDN, DBPDY

REAL(KIND=double)                               :: DBSDT, DBSDN, DBSDY

REAL(KIND=double)                               :: DBUDT, DBUDN, DBUDY

REAL(KIND=double)                               :: DBMUDT, DBMUDN, DBMUDY

REAL(KIND=double)                               :: DMUDT, DMUDN, DMUDY, DSDT, DSDN, DSDY
REAL(KIND=double)                               :: DPDT, DPDN, DPDY, DUDT, DUDN, DUDY

REAL(KIND=double)                               :: NPI, NNI, DNPIDT, DNNIDT

REAL(KIND=double)                               :: DXHDT, DXODT, DXADT

REAL(KIND=double)                               :: VP_IN, VN_IN

REAL(KIND=double)                               :: DUODT, DUIDT, DUADT

REAL(KIND=double)                               :: DBFDT, DBFDY, DBFDN
!
!
!/////////////////////////////////////////////////////////////////////
!
!
!
!
!                   This common block contains all of the thermodynamic
!                   and compositional variables that the EOS code
!                   calculates & it should be included anywhere
!                   these variables need to be referenced
!
!
!     COMMON /EOSVAR/ MUN,MUN_I,MUPROT,MUALFA,MUHAT,
!    1    BU, BS, BPRESS, UTOT, STOT, PTOT,
!    2    XPROT, XNUT, XALFA, XH, X, A,
!    3    BUOUT, BUNUC, BPROUT, BPRNUC, BPRALF, BSOUT,BSNUC,BFOUT,
!    4    BFNUC, BFALFA, FTOT, BFTOT,U_NUC, NSUBI, DNI, RSUBN,
!    5    VSUBN, FBARY, PBARY, MUBARY, MU_MAT,
!    6    XALFA2, XNUT2, XPROT2, IT_NUM
!
!
!                    This common contains the variables that are used
!                    in solving the chemical and pressure equilibrium
!                    equations
!
!     COMMON /EQLEQN/ NNOUT, NPOUT, NOUT, VNOUT, VPOUT, F32_NO, F32_PO,
!    1 MUN_O, MUP_O, ALFDNS, COMPX, SIGMA, SIGSGP, DU, DMU, DUP, DMUP,
!    2 SCRDU, SCRDUP, SCRD, ZETA, FSUBS, FSUBC, MUSUBT, FTRANS, NRATIO,
!    3 VNI, VPI, ETA_NI, ETA_PI, MUP_I, F32_NI, F32_PI, BN, BP, BSUBP,
!    4 H, HPRIM, PSUBI, TSUBC, DHDX, DHDXX, SCRD2, SCRDPP, TFINAL
!
!
!                   This common block contains the derivatives of
!                   the thermodynamic potentials (both baryon and
!                   the totals)
!
!     COMMON /DERIVS/ DBMUDT, DBMUDN, DBMUDY, DMUDT, DMUDN, DMUDY,
!    1                DBPDT, DBPDN, DBPDY, DPDT, DPDN, DPDY,
!    2                DBSDT, DBSDN, DBSDY, DSDT, DSDN, DSDY,
!    3                DBUDT, DBUDN, DBUDY, DUDT, DUDN, DUDY, GAM_S,
!    4                DU_DT, DX_DT, DNI_DT, DEP_DT, DEN_DT,
!    5                DU_DY, DX_DY, DNI_DY, DEP_DY, DEN_DY,
!    6                DU_DN, DX_DN, DNI_DN, DEP_DN, DEN_DN
!
!
!                   This common block contains control flags
!     COMMON /FLAGS/ ADFLAG, NF_FLG, RSFLAG, DBFLAG
!
!
!                   This common block contains the quality of solution
!                   variables that determine how "zeroed" the equilbrium
!                   equations are
!     COMMON /QCHECK/ A1,A2,A3
!
!     double precision ncomp
!     COMMON /DTEST/ DFDOM, DFDL_1,DFDL_2,DFDL_3,DNA_DT,DNA_DN,DNA_DY,
!    1 ncomp
!

END MODULE eos_m4c_module

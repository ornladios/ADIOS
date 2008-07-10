SUBROUTINE eos_nnse_e( j, rho, t, ye, xn, nnc, a_nuc_rep, z_nuc_rep, &
& be_nuc_rep, e_ph, e_elec, e_drip, e_hvy, e_bind, e_no_bind, e_total )
!-----------------------------------------------------------------------
!
!    File:         eos_nnse_e
!    Module:       eos_nnse_e
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/11/04
!
!    Purpose:
!      To compute the various contributions to the internal energy for
!       material not in nse.
!
!    Subprograms called:
!  lectron    : computes electron contribution
!  eos0       : initializes EOS quantities
!  eosnuc_e   : computes nucleon and nuclei contribution
!
!    Input arguments:
!
!  j          : radial zone index
!  rho        : matter density (g/cm3)
!  t          : matter temperature (K)
!  ye         : electron fraction
!  xn(i)      : mass fraction of the ith nucleus
!  nnc        : nuclear abundance array extent
!  a_nuc_rep  : mass number of representative heavy nucleus
!  z_nuc_rep  : charge number of representative heavy nucleus
!  be_nuc_rep : binding energy of representative heavy nucleus (MeV)
!
!    Output arguments:
!  e_ph       : photon energy (ergs g^{-1})
!  e_elec     : electron energy (ergs g^{-1})
!  e_drip     : drip energy (ergs g^{-1})
!  e_hvy      : nuclear internal energy (ergs g^{-1})
!  e_bind     : binding energy (ergs g^{-1})
!  e_no_bind  : total minus binding energy (ergs g^{-1})
!  e_total    : total energy (ergs g^{-1})
!
!    Include files: 
!  kind_module, physcnst.cmn
!  eos_bck_module
!
!----------------------------------------------------------------------c

USE kind_module
USE physcnst_module, ONLY : cm3fm3, rmu, ergmev, kmev, dmnp

USE eos_bck_module, ONLY : dbck, tbck, yebck, ee, ed, eh, b, egy0, jshel

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: nnc           ! composition array dimension

REAL(KIND=double), INTENT(in)    :: rho           ! density (g cm^-3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction
REAL(KIND=double), INTENT(in)    :: a_nuc_rep     ! mass number of representative heavy nucleus
REAL(KIND=double), INTENT(in)    :: z_nuc_rep     ! charge number of representative heavy nucleus
REAL(KIND=double), INTENT(in)    :: be_nuc_rep    ! binding energy of representative heavy nucleus

REAL(KIND=double), INTENT(in), DIMENSION (nnc) :: xn ! nuclaer mass fractions

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: e_ph          ! photon energy (ergs g^{-3})
REAL(KIND=double), INTENT(out)   :: e_elec        ! electron energy (ergs g^{-3})
REAL(KIND=double), INTENT(out)   :: e_drip        ! drip energy (ergs g^{-3})
REAL(KIND=double), INTENT(out)   :: e_hvy         ! nuclear internal energy (ergs g^{-3})
REAL(KIND=double), INTENT(out)   :: e_bind        ! binding energy (ergs g^{-3})
REAL(KIND=double), INTENT(out)   :: e_no_bind     ! total minus binding energy (ergs g^{-3})
REAL(KIND=double), INTENT(out)   :: e_total       ! total energy (ergs g^{-3})

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

REAL(KIND=double), PARAMETER     :: asig = 8.56d-8
REAL(KIND=double)                :: at4           ! photon energy (MeV fm^{-3})
REAL(KIND=double)                :: erad          ! photon energy (MeV nucleon^{-1})

REAL(KIND=double)                :: kfm           ! ( # nucleons/gram )( cm3/fm3 )
REAL(KIND=double)                :: e_cnvrt       ! ( # nucleons/gram )( erg/MeV )

REAL(KIND=double)                :: dbckp         ! density to be passed to eosnuc_e (baryons fm^{-3})
REAL(KIND=double)                :: tbckp         ! temperature to be passed to eosnuc_e (MeV)
REAL(KIND=double)                :: yebckp        ! electron fraction to be passed to eosnuc_e

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize constants
!-----------------------------------------------------------------------

IF ( first ) THEN
  kfm              = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
  e_cnvrt          = ergmev/rmu    ! ( # nucleons/gram )( erg/mev )
  first            = .false.
END IF

!-----------------------------------------------------------------------
!  Convert independent variables
!-----------------------------------------------------------------------

jshel              = j
dbck               = rho * kfm
tbck               = t * kmev
yebck              = ye

!-----------------------------------------------------------------------
!  Photons
!-----------------------------------------------------------------------

at4             = asig * tbck * tbck * tbck * tbck
erad            = at4/dbck

!-----------------------------------------------------------------------
!  Electrons
!-----------------------------------------------------------------------

CALL lectron

!-----------------------------------------------------------------------
!  Nucleons and nuclei
!-----------------------------------------------------------------------

CALL eos0

dbckp           = dbck
tbckp           = tbck
yebckp          = yebck

CALL eosnuc_e( dbckp, tbckp, xn, nnc, a_nuc_rep, z_nuc_rep, be_nuc_rep )

!************************************************ get totals ***********

e_ph            = e_cnvrt * erad
e_elec          = e_cnvrt * ee
e_drip          = e_cnvrt * ed
e_hvy           = e_cnvrt * eh
e_bind          = e_cnvrt * b
e_no_bind       = e_cnvrt * ( eh + erad + ee + ed - egy0 - yebck * dmnp )

e_total         = e_bind + e_no_bind

RETURN
END SUBROUTINE eos_nnse_e

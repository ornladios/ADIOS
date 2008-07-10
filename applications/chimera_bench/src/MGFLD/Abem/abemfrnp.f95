SUBROUTINE abemfrnp( n, e_in, rho, t, xneut, xprot, cmpe, absornp, emitnp )
!-----------------------------------------------------------------------
!
!    File:         abemfrnp
!    Module:       abemfrnp
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute the inverse mean free paths for the absorption
!      and emission of n-type neutrinos on free neutrons and
!      protons. The inverse mean free paths are given in Bruenn,
!      1985, Ap. J. Suupl., 58, 771. Nucleon final state
!      occupancies are not included.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  n           : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  e_in         : neutrino energy (MeV)!       
!  rho         : matter density (g/cm**3)
!  t           : matter temperature (K)
!  xneut       : free neutron mass fraction
!  xprot       : free proton mass fraction
!  cmpe        : electron chemical potential (including rest mass) (MeV)
!
!    Output arguments:
!  absor       : absorption inverse mean free path (/cm)
!  emitnp      : emission inverse mean free path (/cm)!  
!
!    Modules used:
!  kind_module, numerical_module, physcnst_module,
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one
USE physcnst_module, ONLY : cvel, hbar, gv, ga, mp, me, Gw, pi, kmev, dmnp, rmu

USE prb_cntl_module, ONLY : iaefnp, rhoaefnp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)              :: n             ! neutrino flavor index

REAL(KIND=double), INTENT(in)    :: e_in          ! zone centered incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: xneut         ! free neutron mass fraction
REAL(KIND=double), INTENT(in)    :: xprot         ! free proton mass fraction
REAL(KIND=double), INTENT(in)    :: cmpe          ! electron chemical porential

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: absornp       ! inverse mean free path for absorption on free nucleons
REAL(KIND=double), INTENT(out)   :: emitnp        ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

REAL(KIND=double)                :: g2            !  Gw/mp**2 )**2 * ( hbar * cvel )**2
REAL(KIND=double)                :: coef          !  g2 * ( gv**2 + 3.d0 * ga**2 )/pi
REAL(KIND=double)                :: me2           !  electron rest mass squared
REAL(KIND=double)                :: tmev          !  temperature (MeV)
REAL(KIND=double)                :: ron           !  real or effective neutron number
REAL(KIND=double)                :: rop           !  real or effective proton number
REAL(KIND=double)                :: enupdm        !  e_in + dmnp
REAL(KIND=double)                :: enupdm2       !  enupdm^2
REAL(KIND=double)                :: pemitnp       !  ( emission mean free path )/( rop*F_e )
REAL(KIND=double)                :: femitnp       !  ( emission mean free path )/( F_e )
REAL(KIND=double)                :: femitnpl      !  ln(femitnp)
REAL(KIND=double)                :: F_em          !  1 - F_e, 1/( 1 + exp( ( enu + dmnp - cmpe )/t ) )
REAL(KIND=double)                :: F_el          !  ln(F_e)
REAL(KIND=double)                :: eta           !  ( e_in + dmnp - cmpe )/tmev
REAL(KIND=double)                :: ex            !  exp(eta)
REAL(KIND=double)                :: cmppos        !  positron chemical potential
REAL(KIND=double)                :: F_p           !  1/( 1 + exp( ( e_in - dmnp - cmppos )/t ) )
REAL(KIND=double)                :: enumdm        !  e_in - dmnp
REAL(KIND=double)                :: enumdm2       !  enumdm^2

REAL(KIND=double)                :: emitnpl       !  ln(emitnp)

REAL(KIND=double), PARAMETER     :: expmax = 3.d02

REAL(KIND=double)                :: fexp

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  .Return if iaefnp = 0 or rho  >  rhoaefnp or n > 2
!-----------------------------------------------------------------------

IF ( iaefnp == 0  .or.  rho  >  rhoaefnp  .or.  n > 2 ) THEN
  emitnp             = zero
  absornp            = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Set emitnp and absornp to zero and return if both xneut and
!   xprot are zero.
!-----------------------------------------------------------------------

IF ( xneut == zero  .and.  xprot == zero ) THEN
  emitnp             = zero
  absornp            = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Initialize.
!
!  g2   : in units of cm**2/MeV**2.
!  coef : coefficient of the emission and absornpption inverse mean free paths.
!-----------------------------------------------------------------------

IF ( first ) THEN
  first              = .false.
  g2                 = ( Gw/mp**2 )**2 * ( hbar * cvel )**2
  coef               = g2 * ( gv**2 + 3.d0 * ga**2 )/pi
  me2                = me * me
END IF

!-----------------------------------------------------------------------
!  ron : number of free neutrons per unit volume (/cm**3).
!  rop : number of free protons per unit volume (/cm**3).
!-----------------------------------------------------------------------

tmev                 = kmev * t
ron                  = xneut * rho/rmu
rop                  = xprot * rho/rmu

!-----------------------------------------------------------------------
!
!    \\\\\ ABSORPTION AND EMISSION INV MFP'S FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  pemitnp: ( emission mean free path )/( rop*F_e ).
!  femitnp: ( emission mean free path )/( F_e ).
!-----------------------------------------------------------------------

IF ( n == 1 ) THEN

  enupdm             = e_in + dmnp
  enupdm2            = enupdm * enupdm
  pemitnp            = coef * enupdm2 * dsqrt( one - me2/enupdm2 )
  femitnp            = rop * pemitnp

!-----------------------------------------------------------------------
!  Case that the free proton fraction is zero
!-----------------------------------------------------------------------

  IF ( femitnp == zero ) THEN
    emitnp           = zero
    eta              = ( e_in + dmnp - cmpe )/tmev
    ex               = fexp(eta)
    F_em             = ex/( one + ex )
    absornp          = ron * pemitnp * F_em
    RETURN
  END IF

!-----------------------------------------------------------------------
!  General case compute logs of inv mfp's
!-----------------------------------------------------------------------

  femitnpl           = dlog(femitnp)
  eta                = ( e_in + dmnp - cmpe )/tmev
  IF ( eta > expmax ) THEN
    F_el             = -eta
  ELSE
    ex               = fexp(eta) 
    F_el             = -DLOG( one + ex )
  END IF

!-----------------------------------------------------------------------
!  Log of emission inverse mean free path; emission inverse mean free
!   path.
!-----------------------------------------------------------------------

  emitnpl            = femitnpl + F_el
  emitnp             = fexp(emitnpl)

!-----------------------------------------------------------------------
!  Compute the asorption inverse mean free path
!-----------------------------------------------------------------------

  F_em               = ex/( one + ex )
  absornp            = ron * pemitnp * F_em

  RETURN
  
END IF

!-----------------------------------------------------------------------
!
!  \\\\\ ABSORPTION AND EMISSION INV MFP'S FOR E-ANTINEUTRINOS /////
!
!-----------------------------------------------------------------------

IF ( n == 2 ) THEN

  IF ( e_in - me <= DABS(dmnp) ) THEN
    emitnp           = zero
    absornp          = zero
    RETURN
  END IF

!-----------------------------------------------------------------------
!  Absorption and emission inv mfp's for e-antineutrinos
!
!  pemitnp: ( emission mean free path )/( ron*F_p )
!-----------------------------------------------------------------------

  enumdm             = e_in - dmnp
  enumdm2            = enumdm * enumdm
  pemitnp            = coef * enumdm2 * dsqrt( one - me2/enumdm2 )

!-----------------------------------------------------------------------
!  Positron occupation number for positrons of energy e_in - dmnp
!-----------------------------------------------------------------------

  cmppos             = -cmpe
  eta                = ( e_in - dmnp - cmppos )/tmev
  ex                 = fexp(eta)
  F_p                = one/( one + ex )

!-----------------------------------------------------------------------
!  inv mfp's
!-----------------------------------------------------------------------

  emitnp             = ron * pemitnp * F_p
  absornp            = rop * pemitnp * ( one - F_p )

  RETURN
  
END IF

STOP
END SUBROUTINE abemfrnp

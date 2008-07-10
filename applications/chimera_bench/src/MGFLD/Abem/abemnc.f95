SUBROUTINE abemnc( n, e_in, rho, t, xh, ah, zh, cmpn, cmpp, cmpe, absrnc, emitnc )
!-----------------------------------------------------------------------
!
!    File:         abemnc
!    Module:       abemnc
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To computes the inverse mean free path for the absorption of e-type neutrinos
!       on nuclei and the inverse process (emission of e-type neutrinos
!       by electron capture on nuclei). The Fuller, Fowler, Neuman 1982, Ap. J. 252, 715
!       approximation for the heavy nucleus matrix element as given in Bruenn 1985,
!       Ap. J. Suupl., 58, 771 is used.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  n           : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  e_in        : neutrino energy (MeV)
!  rho         : matter density (g/cm**3)
!  t           : matter temperature (K)
!  xh          : heavy nuclei mass fraction
!  ah          : heavy nuclei mass number
!  zh          : heavy nuclei charge number
!  cmpn        : free neutron chemical potential (excluding rest mass) (MeV)
!  cmpp        : free proton chemical potential (excluding rest mass) (MeV)
!  cmpe        : electron chemical potential (including rest mass) (MeV)
!                 between the exitation energy of daughter and parent nucleus (MeV)
!
!    Output arguments:
!  absrnc      : absorption inverse mean free path (/cm)
!  emitnc      : emission inverse mean free path (/cm)
!
!    Modules used:
!  kind_module, numerical_module, physcnst_module,
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, one, epsilon
USE physcnst_module, ONLY : Gw, mp, hbar, cvel, pi, ga, kmev, me, dmnp, rmu

USE prb_cntl_module, ONLY : iaence, edmpe

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)              :: n             ! neutrino flavor index

REAL(KIND=double), INTENT(in)    :: e_in          ! zone centered incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: xh            ! heavy nuclei mass fraction
REAL(KIND=double), INTENT(in)    :: ah            ! heavy nuclei mass number
REAL(KIND=double), INTENT(in)    :: zh            ! heavy nuclei charge number
REAL(KIND=double), INTENT(in)    :: cmpn          ! neutron chemical porential
REAL(KIND=double), INTENT(in)    :: cmpp          ! proton chemical porential
REAL(KIND=double), INTENT(in)    :: cmpe          ! electron chemical porential

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout) :: absrnc        ! inverse mean free path for absorption on free nucleons
REAL(KIND=double), INTENT(inout) :: emitnc        ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

REAL(KIND=double)                :: g2            ! square of the Fermi constant in units of cm^{2} MeV^{-2}
REAL(KIND=double)                :: coef          ! coefficient of the emis and abs inverse mean free paths (cm^{-1})
REAL(KIND=double)                :: me2           ! square of the electron rest mass
REAL(KIND=double)                :: nh            ! number of neutrons
REAL(KIND=double)                :: coefn         ! number of neutron holes
REAL(KIND=double)                :: zm20          ! number of available protons
REAL(KIND=double)                :: eelec         ! electron energy (MeV)
REAL(KIND=double)                :: eelec2        ! electron energy squared (MeV)
REAL(KIND=double)                :: ronc          ! number density of nuclei (cm^{-3})
REAL(KIND=double)                :: ceelec        ! coefficient for computing the emission inverse mean free path
REAL(KIND=double)                :: pemit         ! emission inverse mean free path per electron
REAL(KIND=double)                :: eta           ! (eelec - cmpe)/tmev
REAL(KIND=double)                :: F_e           ! electron occupation number
REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double)                :: etam          ! ( e_in + dmnp + cmpn - cmpp - cmpe )/tmev
REAL(KIND=double)                :: fexp          ! exponential function

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set rates to zero and return if iaence = 0 or if n ne 1
!-----------------------------------------------------------------------

IF ( iaence == 0  .or.  n /= 1 ) THEN
  emitnc           = zero
  absrnc           = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.
  g2               = ( Gw/mp**2 )**2 * ( hbar * cvel )**2
  coef             = ( 2.d0/7.d0 ) * ( g2/pi ) * ( ga**2 )
  me2              = me * me
END IF

!-----------------------------------------------------------------------
!  Compute coefn, the number of neutron holes
!-----------------------------------------------------------------------

nh                 = ah - zh
coefn              = DMIN1( 40.d+00 - nh, 6.d+00 )
IF ( coefn < zero ) THEN
  emitnc           = zero
  absrnc           = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  The Bowers and Wilson prescription is as follows:
!
!  if (coefn .lt. 6./16.) coefn = 6./16.
!  coefn        = coefn/5.
!
!  Compute the number of available protons.
!-----------------------------------------------------------------------

zm20               = DMIN1( zh - 20.d+00, 8.d+00 )
IF ( zm20 <= zero ) THEN
  emitnc           = zero
  absrnc           = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Compute the electron energy necessary to produce an e-neutrino
!   of energy e_in.
!-----------------------------------------------------------------------

eelec              = e_in + cmpn - cmpp + edmpe + dmnp
eelec2             = eelec * eelec

!-----------------------------------------------------------------------
!  Compute the emission rate, pemit, per electron
!-----------------------------------------------------------------------

ronc               = xh * rho/( ah * rmu )
ceelec             = DSQRT( DMAX1( one - me2/eelec2, epsilon ) )
pemit              = coef * ronc * zm20 * eelec2 * ceelec

!-----------------------------------------------------------------------
!  Compute the electron occupation number, F_e
!-----------------------------------------------------------------------

tmev               = kmev * t
eta                = ( eelec - cmpe )/tmev
F_e                = one/( one + fexp(eta) )

!-----------------------------------------------------------------------
!  Compute the emission inverse mean free path
!-----------------------------------------------------------------------

emitnc             = pemit * F_e * coefn

!-----------------------------------------------------------------------
!  Compute the absorption inverse mean free path using detailed
!   balance.
!-----------------------------------------------------------------------

etam               = ( e_in + dmnp + cmpn - cmpp - cmpe )/tmev
absrnc             = emitnc * fexp(etam)


RETURN
END SUBROUTINE abemnc

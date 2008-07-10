SUBROUTINE abemhnc( n, ie, rho, t, xh, ah, cmpn, cmpp, cmpe, absrnc, &
& emitnc )
!-----------------------------------------------------------------------
!
!    File:         abemhnc
!    Module:       abemhnc
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute the inverse mean free paths for the absorption and emission
!       of n-type neutrinos on certain nuclei using cross sections provided
!       by Wick Haxton.
!
!    Variables that must be passed through common:
!
!        rncem(k,kp): e-neutrino - nucleus absorption rate (per iron-like nucleus) and per
!          incident neutrino of energy k and final electrons of energy kp as given by Haxton.
!         
!
!        rncep(k,kp): e-antineutrino - nucleus absorption rate (per iron-like nucleus) and per
!          incident neutrino of energy k and final electrons of energy kp as given by Haxton.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  n           : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  ie          : neutrino energy zone index   
!  rho         : matter density (g/cm**3)
!  t           : matter temperature (K)
!  xh          : heavy nucleus mass fraction
!  ah          : heavy nucleus mass number
!  cmpn        : free neutron chemical potential (excluding rest mass) (MeV)
!  cmpp        : free proton chemical potential (excluding rest mass) (MeV)
!  cmpe        : electron chemical potential (including rest mass) (MeV)
!
!    Output arguments:
!  absrnc      :  absorption inverse mean free path (/cm)
!  emitnc      :  emission inverse mean free path (/cm)
!
!    Input arguments (common):
!
!  iaencnu     : 0; inverse mean free paths set to zero.
!                1; inverse mean free paths computed.
!  roaencnu    : density above which rates are set to zero.
!
!    Modules used:
!  kind_module, numerical_module, physcnst_module,
!  abem_module, nu_energy_grid_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, one
USE physcnst_module, ONLY : rmu, kmev, dmnp

USE abem_module, ONLY : rncem, rncep
USE nu_energy_grid_module, ONLY : nnugp, unui
USE prb_cntl_module, ONLY : iaencnu, roaencnu

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: n             ! neutrino flavor index
INTEGER, INTENT(in)              :: ie            ! neutrino energy index

REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: xh            ! heavy nuclei mass fraction
REAL(KIND=double), INTENT(in)    :: ah            ! heavy nuclei mass number
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

INTEGER                          :: je            ! do index

REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double)                :: xnuc          ! number density of nuclei (cm^{-3})
REAL(KIND=double)                :: eelec         ! electron energy (MeV)
REAL(KIND=double)                :: eta           ! ( eelec - cmpe )/tmev
REAL(KIND=double)                :: ex            ! exp(eta)
REAL(KIND=double)                :: F_em          ! electron blocking factor
REAL(KIND=double)                :: absrnct       ! term in the absorption inverse mean free path
REAL(KIND=double)                :: eptron        ! positron energy (MeV)
REAL(KIND=double)                :: cmppos        ! positron chemical potential (MeV)
REAL(KIND=double)                :: etam          ! ( cmpe - cmpn + cmpp - dmnp - unui(ie) )/tmev
REAL(KIND=double)                :: F_pm          ! positron blocking factor
REAL(KIND=double)                :: fexp          ! exponential function

EXTERNAL fexp

!-----------------------------------------------------------------------
!  Set emitnc and absrnc to zero and return if
!               iaencnu = 0
!  or
!               rho > roaencnu
!  or
!               n = 3
!-----------------------------------------------------------------------

IF ( iaencnu == 0  .or.  rho  >  roaencnu  .or.  n == 3 ) THEN
  emitnc           = zero
  absrnc           = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

emitnc            = zero
absrnc            = zero
tmev              = kmev * t
xnuc              = ( xh/( rmu * ah ) ) * rho

!-----------------------------------------------------------------------
!
!       \\\\\ E-NEUTRINO-NUCLEUS ABSORPTION AND EMISSION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  The neutrino has energy unu(ie); a sum is made over the possible
!   final electron enrgies unui(je).
!-----------------------------------------------------------------------

IF ( n == 1 ) THEN

  DO je = 1,nnugp(1)
    IF ( rncem(ie,je) /= zero ) THEN

!-----------------------------------------------------------------------
!  Electron energy
!-----------------------------------------------------------------------

      eelec        = unui(je)

!-----------------------------------------------------------------------
!  Compute 1 -  electron occupation number) = F_em
!-----------------------------------------------------------------------

      eta          = ( eelec - cmpe )/tmev
      ex           = fexp(eta)
      F_em         = ex/( one + ex )

!-----------------------------------------------------------------------
!  Compute the absorption rate
!-----------------------------------------------------------------------

      absrnct      = xnuc * F_em * rncem(ie,je)
      absrnc       = absrnc + absrnct

!-----------------------------------------------------------------------
!  Compute the emission rate using detailed balance
!-----------------------------------------------------------------------

      etam         = ( cmpe - cmpn + cmpp - dmnp - unui(ie) )/tmev
      emitnc       = emitnc + absrnct * fexp(etam)

    END IF ! ( rncem(ie,je) /= 0 )

  END DO

END IF ! ( n == 1 )

!-----------------------------------------------------------------------
!
!     \\\\\ E-ANTINEUTRINO-NUCLEUS ABSORPTION AND EMISSION /////
!
!-----------------------------------------------------------------------

IF ( n == 2 ) THEN

  DO je = 1,nnugp(2)

    IF ( rncep(ie,je) /= 0 ) THEN

!-----------------------------------------------------------------------
!  Positron energy
!-----------------------------------------------------------------------

      eptron       = unui(je)

!-----------------------------------------------------------------------
!  Compute (1 -  positron occupation number) = F_em
!-----------------------------------------------------------------------

      cmppos       = -cmpe
      eta          = ( eptron - cmppos )/tmev
      ex           = fexp(eta)
      F_pm         = ex/( one + ex )

!-----------------------------------------------------------------------
!  Compute the absorption rate
!-----------------------------------------------------------------------

      absrnct      = xnuc * F_pm * rncep(ie,je)
      absrnc       = absrnc + absrnct

!-----------------------------------------------------------------------
!  Compute the emission rate using detailed balance
!-----------------------------------------------------------------------

      etam         = ( cmppos + cmpn - cmpp + dmnp - unui(ie) )/tmev
      emitnc       = emitnc + absrnct*fexp(etam)

    END IF ! ( rncep(ie,je) /= 0 )

  END DO

END IF ! ( n == 2 )


RETURN
END SUBROUTINE abemhnc

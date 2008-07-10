SUBROUTINE e_degen( rho, t, ye, e_eta )
!-----------------------------------------------------------------------
!
!    File:         e_degen
!    Module:       e_degen
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/13/01
!
!    Purpose:
!      To compute the degeneracy parameter of the electrons
!
!    Variables that must be passed through common:
!
!
!    Subprograms called:
!      tadv0nr, tadv1nr, tadv0, tadv1, tadv4
!
!    Input arguments:
!  rho         : matter density of zone j (g/cm**3)
!  t           : matter temperature of zone j (g/cm**3)
!  ye          : matter electron fraction
!
!    Output arguments:
!  lconv: true : iteration for the temperature, density, and
!                 relativistic parameters has converged
!        false : iteration for the temperature, density, and
!                 relativistic parameters has not converged
!
!    Input arguments (common):
!  rho(j)      : matter density of zone j (g/cm**3)
!  t(j)        : matter temperature of zone j (g/cm**3)
!  ye(j)       : matter electron fraction
!
!    Output arguments (common):
!        none
!
!    Include files:
!      kind_module, numerical_module, physcnst_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : third
USE physcnst_module, ONLY : pi, cm3fm3, rmu, hbarc, me, kmev

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)    :: rho           ! density (g cm^{-3})
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: e_eta         ! (electron Fermi energy)/KT

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

REAL(KIND=double)                :: pi2           ! pi^{2}
REAL(KIND=double)                :: kfm           ! converts g/cm**3 to number/fm**3
REAL(KIND=double)                :: rhofm         ! density (number per fm^{3})
REAL(KIND=double)                :: pfc           ! electron Fermi momentum
REAL(KIND=double)                :: efk           ! electron Fermi energy

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Initialize constants...........................................

IF ( first ) THEN
  pi2              = pi**2
  kfm              = cm3fm3/rmu
  first            = .false.
END IF

!........Electron degeneracy parameter, etatst.........................

rhofm              = rho * kfm
pfc                = hbarc * ( 3.d0 * pi2 * rhofm * ye )**third

IF ( rhofm * ye > 1.d-10 ) THEN
  efk              = DSQRT(pfc * pfc + me * me) - me
ELSE
  efk              = pfc * pfc/( 2.d0 * me )
END IF ! rhofm * ye(j) > 1.d-10

e_eta              = efk/( t * kmev )

RETURN
END SUBROUTINE e_degen

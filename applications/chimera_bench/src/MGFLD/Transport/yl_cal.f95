SUBROUTINE yl_cal( jr_min, jr_max, rho, ye, psi0, yenu, yenubar, yxnu, &
& yxnubar, yl, nx, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         yl_cal
!    Module:       yl_cal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/00
!
!    Purpose:
!      To compute the lepton fraction.
!
!    Variables that must be passed through common:
!  agr_nu(j) : GR lapse function
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min    : inner zone for which the calculation of yl is to be made
!  jr_max    : outer zone for which the calculation of yl is to be made
!  rho       : density (g xm^{-3})
!  ye        : electron fraction
!  psi0      : zeroth angular moment of the NDS
!  nx        : x_array extent
!  nez       : neutrino energy array extent
!  nnu       : neutrino flavor array extent
!
!    Output arguments:
!  yl        : lepton fraction
!  yenu      : electron neutrino fraction
!  yenubar   : electron antineutrino fraction
!  yxnu      : muon or tau neutrino fraction
!  yxnubar   : muon or tau antineutrino fraction
!
!    Input arguments (common):
!  nnugp(n)  : number of neutrino energy groups for n-neutrinos
!  unu(k)    : zone centered energy of neutrino energy zone k
!  dunu(i)   : unu(j,k+1) - unu(j,k)
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, ncoef
USE physcnst_module, ONLY : rmu

USE nu_dist_module, ONLY : unu, dunu
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)               :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)               :: nx            ! x-array extent
INTEGER, INTENT(in)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: rho         ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: ye          ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu) :: psi0        ! zeroth angular moment of the NDS

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx)        :: yl          ! lepton fraction
REAL(KIND=double), INTENT(out), DIMENSION(nx)        :: yenu        ! electron neutrino fraction
REAL(KIND=double), INTENT(out), DIMENSION(nx)        :: yenubar     ! electron antineutrino fraction
REAL(KIND=double), INTENT(out), DIMENSION(nx)        :: yxnu        ! muon or tau neutrino fraction
REAL(KIND=double), INTENT(out), DIMENSION(nx)        :: yxnubar     ! muon or tau antineutrino fraction

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index

REAL(KIND=double), DIMENSION(2)   :: n_nu          ! neutrino number density

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                   \\\\\ LEPTON FRACTION /////
!
!-----------------------------------------------------------------------

DO j = jr_min,jr_max

  DO n = 1,2
    n_nu(n)        = zero
    IF ( nnugp(n) /= 0 ) THEN
      DO k = 1,nnugp(n)
        n_nu(n)    = n_nu(n) + ( ncoef/rho(j) ) * unu(j,k) * unu(j,k) * dunu(j,k) * psi0(j,k,n)
      END DO
    END IF !  nnugp(n) ne 0
  END DO

  yenu(j)          = n_nu(1) * rmu
  yenubar(j)       = n_nu(2) * rmu
  yxnu(j)          = zero
  yxnubar(j)       = zero
  yl(j)            = ye(j) + yenu(j) - yenubar(j)

END DO

RETURN
END SUBROUTINE yl_cal

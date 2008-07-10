SUBROUTINE nu_energy_advct_z( kmin, kmax, ki_ray, kj_ray, i_radial)
!-----------------------------------------------------------------------
!
!    File:         nu_energy_advct_z
!    Module:       nu_energy_advct_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/07
!
!    Purpose:
!      To advance the advance the radii, densities, temperatures,
!       gravitational masses, and GR gammas along with the neutrino
!       distributions due to energy advection.
!
!    Variables that must be passed through common:
!
!    Subprograms called:
!  e_advct_z        : performs the neutrino energy advection due to z-hydro
!  rebal            : enforces the exclusion principle
!
!    Input arguments:
!  kmin             : inner z-zone for which energy advection is to be calculated
!  kmax             : outer z-zone for which energy advection is to be calculated
!  ki_ray           : x (radial) index of a specific z (azimhthal) ray
!  kj_ray           : y (angular) index of a specific z (azimhthal) ray
!  i_radial         : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!
!    Output arguments:
!        none
!
!    Output arguments (common):
!  psi0             : zero angular moments of the neutrino occupation number
!  unujv(j,n,j_ray) : energy transferred to n-type neutrinos by energy advection (ergs)
!  unucrv(n,j_ray)  : total energy transferred to n-type neutrinos by energy advection (ergs)
!  dunujvdt(j,n)    : rate of energy transferred to n-type neutrinos by energy advection
!  dndt_v(j,k,n)    : net rate of n-neutrino production by energy advection
!
!    Include files:
!  kind_module, array_module, numerical_module
!  e_advct_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nnu
USE numerical_module, ONLY : zero

USE e_advct_module, ONLY : psi0, psi0_a, dpsivmx, psivmin, &
& jdt_nu_e_advct
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)          :: kmin      ! minimum z-zone index
INTEGER, INTENT(in)          :: kmax      ! maximum z-zone index
INTEGER, INTENT(in)          :: ki_ray    ! x (radial) index of a specific z (azimhthal) ray
INTEGER, INTENT(in)          :: kj_ray    ! y (angular) index of a specific z (azimhthal) ray
INTEGER, INTENT(in)          :: i_radial  ! the unshifted radial zone corresponding to ki_ray, kj_ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                      :: j         ! radial zone index
INTEGER                      :: k         ! neutrino energy index
INTEGER                      :: n         ! neutrino flavor index

REAL(KIND=double)            :: denu      ! maximum relative change of psi0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ADVECT THE NEUTRINOS IN ENERGY /////
!
!  Advance neutrino occupation probabilities due to the energy
!   advection of neutrinos.
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  CALL e_advct_z( kmin, kmax, n, ki_ray, kj_ray, i_radial )
END DO

!-----------------------------------------------------------------------
!  Restore updated neutrino distribution to psi0 array and determine
!   the largest relative change of psi0.
!-----------------------------------------------------------------------

dpsivmx             = zero
DO j = kmin,kmax
  DO k = 1,nnugpmx
    DO n = 1,nnu
    IF ( nnugp(n) == 0 ) CYCLE
      denu          = DABS( psi0_a(j,k,n) - psi0(j,k,n) )/( psi0_a(j,k,n) + psivmin(n) )
      IF ( denu > dpsivmx(n) ) THEN
        dpsivmx(n)  = denu
        jdt_nu_e_advct(n) = j
      END IF
      psi0(j,k,n)   = psi0_a(j,k,n)
    END DO !  n = 1,nnu
  END DO !  k = 1,nnugpmx
END DO !  j = kmin,kmax

!-----------------------------------------------------------------------
!  Prevent overfilling of neutrinos states
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  CALL rebal( kmin, kmax, n )
END DO

RETURN
END SUBROUTINE nu_energy_advct_z
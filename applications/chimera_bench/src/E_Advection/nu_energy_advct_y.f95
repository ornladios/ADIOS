SUBROUTINE nu_energy_advct_y( jmin, jmax, ji_ray, jk_ray, i_radial)
!-----------------------------------------------------------------------
!
!    File:         nu_energy_advct_y
!    Module:       nu_energy_advct_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/16/01
!
!    Purpose:
!      To advance the advance the radii, densities, temperatures,
!       gravitational masses, and GR gammas along with the neutrino
!       distributions due to energy advection.
!
!    Variables that must be passed through common:
!
!    Subprograms called:
!  e_advct_y        : performs the neutrino energy advection due to y-hydro
!  rebal            : enforces the exclusion principle
!
!    Input arguments:
!  jmin             : inner y-zone for which energy advection is to be calculated
!  jmax             : outer y-zone for which energy advection is to be calculated
!  ji_ray           : x (radial) index of a specific y (angular) ray
!  jk_ray           : z (azimuthal) index of a specific y (angular) ray
!  i_radial         : the unshifted radial zone (angular ray) corresponding to ji_ray, jk_ray
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

INTEGER, INTENT(in)          :: jmin      ! minimum y-zone index
INTEGER, INTENT(in)          :: jmax      ! maximum y-zone index
INTEGER, INTENT(in)          :: ji_ray    ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)          :: jk_ray    ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)          :: i_radial  ! the unshifted radial zone corresponding to ji_ray, jk_ray

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
  CALL e_advct_y( jmin, jmax, n, ji_ray, jk_ray, i_radial )
END DO

!-----------------------------------------------------------------------
!  Restore updated neutrino distribution to psi0 array and determine
!   the largest relative change of psi0.
!-----------------------------------------------------------------------

dpsivmx             = zero
DO j = jmin,jmax
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
END DO !  j = jmin,jmax

!-----------------------------------------------------------------------
!  Prevent overfilling of neutrinos states
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  CALL rebal( jmin, jmax, n )
END DO

RETURN
END SUBROUTINE nu_energy_advct_y
SUBROUTINE nu_energy_advct_y( jmin, jmax, j_ray, j_ray_dim )
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
!  eddington        : computes the neutrino eddington factor for a given energy group
!  e_advct          : performs the neutrino energy advection
!  rebal            : enforces the exclusion principle
!
!    Input arguments:
!  jmin             : inner zone for which energy advection is to be calculated
!  jmax             : outer zone for which energy advection is to be calculated
!  j_ray            : index denoting a specific angular ray
!  j_ray_dim        : number of angular rays to put on a processor
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
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)          :: jmin      ! minimum zone index
INTEGER, INTENT(in)          :: jmax      ! maximum zone index
INTEGER, INTENT(in)          :: j_ray     ! index denoting a specific angular ray
INTEGER, INTENT(in)          :: j_ray_dim ! number of angular rays assigned to a processor

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
!        Advance neutrino occupation probabilities due to the
!         energy advection of neutrinos.
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  print *,' nu_energy_advct_y 1, n=',n
  CALL e_advct_y( jmin, jmax, n, j_ray, j_ray_dim )
END DO

!-----------------------------------------------------------------------
!        Restore updated neutrino distribution to psi0 array and
!         determine the largest relative change of psi0.
!-----------------------------------------------------------------------

DO n = 1,nnu
  print *,' nu_energy_advct_y 2, n=',n
  IF ( nnugp(n) == 0 ) CYCLE
  dpsivmx(n)        = zero
  DO k = 1,nnugp(n)
    DO j = jmin,jmax
      denu          = DABS( psi0_a(j,k,n) - psi0(j,k,n) )/( psi0_a(j,k,n) + psivmin(n) )
      IF ( denu > dpsivmx(n) ) THEN
        dpsivmx(n)  = denu
        jdt_nu_e_advct(n) = j
      END IF
      psi0(j,k,n)   = psi0_a(j,k,n)
    END DO
  END DO
END DO

!........Prevent overfilling of neutrinos states........................

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  print *,' nu_energy_advct_y 3, n=',n
  CALL rebal( jmin, jmax, n )
END DO

RETURN
END SUBROUTINE nu_energy_advct_y
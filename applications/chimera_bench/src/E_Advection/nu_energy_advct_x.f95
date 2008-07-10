SUBROUTINE nu_energy_advct_x( jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         nu_energy_advct_x
!    Module:       nu_energy_advct_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/16/01
!
!    Purpose:
!      To advance the advance the the neutrino distributions due to neutrino 
!       energy advection following an x-hydro step
!
!    Subprograms called:
!  eddington                : computes the neutrino eddington factor for a given energy group
!  e_advct                  : performs the neutrino energy advection
!  rebal                    : enforces the exclusion principle
!
!    Input arguments:
!  jr_min                   : inner zone for which energy advection is to be calculated
!  jr_max                   : outer zone for which energy advection is to be calculated
!  ij_ray                   : j-index of a radial ray
!  ik_ray                   : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Output arguments (common):
!  psi0                     : zero angular moments of the neutrino occupation number
!  unujv(j,n,ij_ray,ik_ray) : energy transferred to n-type neutrinos by energy advection (ergs)
!  unucrv(n,ij_ray,ik_ray)  : total energy transferred to n-type neutrinos by energy advection (ergs)
!  dunujvdt(j,n)            : rate of energy transferred to n-type neutrinos by energy advection
!  dndt_v(j,k,n)            : net rate of n-neutrino production by energy advection
!
!    Include files:
!  kind_module, array_module, numerical_module
!  e_advct_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nnu
USE numerical_module, ONLY : zero

USE e_advct_module, ONLY : r, ra, psi0, psi0_a, dpsivmx, psivmin, &
& jdt_nu_e_advct
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)      :: jr_min        ! minimum zone index
INTEGER, INTENT(in)      :: jr_max        ! maximum zone index
INTEGER, INTENT(in)      :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)      :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                  :: j             ! radial zone index
INTEGER                  :: k             ! neutrino energy index
INTEGER                  :: n             ! neutrino flavor index

REAL(KIND=double)        :: denu          ! maximum relative change of psi0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute the flux and Eddington factors
!-----------------------------------------------------------------------

CALL eddington( jr_min, jr_max )

!-----------------------------------------------------------------------
!
!             \\\\\ ADVECT THE NEUTRINOS IN ENERGY /////
!
!  Advance neutrino occupation probabilities due to the energy
!   advection of neutrinos.
!-----------------------------------------------------------------------

ra(jr_max+1)          = ra(jr_max) + ( ra(jr_max) - ra(jr_max-1) )
r (jr_max+1)          = r (jr_max) + ( r (jr_max) - r (jr_max-1) )

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  CALL e_advct_x( jr_min, jr_max, n, ij_ray, ik_ray )
END DO

!-----------------------------------------------------------------------
!  Restore updated neutrino distribution to psi0 array and determine
!   determine the largest relative change of psi0.
!-----------------------------------------------------------------------

dpsivmx             = zero
DO j = jr_min,jr_max
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
END DO !  j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Prevent overfilling of neutrino states
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  CALL rebal( jr_min, jr_max, n )
END DO

RETURN
END SUBROUTINE nu_energy_advct_x
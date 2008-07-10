SUBROUTINE nu_energy_agr_advct_x( jr_min, jr_max )
!-----------------------------------------------------------------------
!
!    File:         nu_energy_agr_advct_x
!    Module:       nu_energy_agr_advct_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/19/06
!
!    Purpose:
!      To advance the advance the the neutrino distributions due to
!       neutrino energy advection following an x-hydro step
!
!    Subprograms called:
!  e_advct_agr_x : performs the neutrino energy advection
!  rebal         : enforces the exclusion principle
!
!    Input arguments:
!  jr_min        : inner zone for which energy advection is to be calculated
!  jr_max        : outer zone for which energy advection is to be calculated
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  e_advct_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nnu
USE numerical_module, ONLY : zero, ncoef, ecoef

USE e_advct_module, ONLY : agrajmh, ncoefaa, ecoefaa
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx,  unui, dunui, unubi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)      :: jr_min        ! minimum zone index
INTEGER, INTENT(in)      :: jr_max        ! maximum zone index

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                  :: j             ! shifted radial zone index
INTEGER                  :: n             ! neutrino flavor index

REAL(KIND=double)        :: agrajmh_3     ! 1/agrajmh^{3}
REAL(KIND=double)        :: agrajmh_4     ! 1/agrajmh^{4}

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
  CALL e_advct_agr_x( jr_min, jr_max, n )
END DO ! n = 1,nnu

!-----------------------------------------------------------------------
!  Compute neutrino number and energy coefficients
!-----------------------------------------------------------------------

DO j = jr_min, jr_max
  agrajmh_3             = 1.d0/agrajmh(j)**3
  agrajmh_4             = 1.d0/agrajmh(j)**4
  ncoefaa(j,1:nnugpmx)  = ncoef * unui(1:nnugpmx)**2 * dunui(1:nnugpmx) * agrajmh_3
  ecoefaa(j,1:nnugpmx)  = ecoef * unui(1:nnugpmx)**3 * dunui(1:nnugpmx) * agrajmh_4
END DO ! j = jr_min, jr_max

!-----------------------------------------------------------------------
!  Prevent overfilling of neutrinos states
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  CALL rebal( jr_min, jr_max, n )
END DO

RETURN
END SUBROUTINE nu_energy_agr_advct_x
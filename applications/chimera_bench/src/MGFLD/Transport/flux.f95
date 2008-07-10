SUBROUTINE flux( jr_min, jr_max, n )
!-----------------------------------------------------------------------
!
!    File:         flux
!    Module:       flux
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/22/04
!
!    Purpose:
!      To calculate at neutrino flux.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  jr_min         : inner zone for which calculation w is to be made
!  jr_max         : outer zone for which calculation w is to be made
!  n              : neutrino flavor index
!
!    Output arguments:
!      none
!
!    Input arguments (common):
!  nnugp(n)       : number of neutrino energy groups for n-neutrinos
!  ecoefae(j,k)   : energy sum neutrino states per unit volume at radial zone j energy unu(j,k)
!  psi1(j,k,n)    : first angular moment of the neutrino distribution function
!
!    Output arguments (common):
!  fluxnuk(j,k,n) :  rate of energy flow through radial zone j due to
!                     neutrinos of type n in group k [ergs sec^{-1} cm^{-2}].
!  fluxnu(j,n)    :  rate of energy flow through radial zone j due to
!                     neutrinos of type n [ergs sec^{-1} cm^{-2}].
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, third, half
USE physcnst_module, ONLY : cvel

USE nu_dist_module, ONLY : ecoefae, stwt, psi1, fluxnuk, fluxnu
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: n             ! neutrino flavor index

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

fluxnu(jr_min:jr_max,n)             = zero
IF ( nnugp(n) == 0 ) RETURN
fluxnuk(jr_min:jr_max,1:nnugp(n),n) = cvel * psi1(jr_min:jr_max,1:nnugp(n),n) &
&                                   * ecoefae(jr_min:jr_max,1:nnugp(n)) * stwt(n)

DO j = jr_min,jr_max
  fluxnu(j,n)                       = SUM( fluxnuk(j,1:nnugp(n),n) )
END DO ! j = jr_min,jr_max

RETURN
END SUBROUTINE flux

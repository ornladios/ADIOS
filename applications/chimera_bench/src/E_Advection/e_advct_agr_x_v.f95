SUBROUTINE e_advct_agr_x_v( j, n )
!-----------------------------------------------------------------------
!
!    File:         e_advct_agr_x_v
!    Module:       e_advct_agr_x_v
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/19/06
!
!    Purpose:
!      To compute the neutrino advection velocity through the
!       energy grid.
!
!    Subprograms called:
!  eddington     : computes the ration of psi2 to psi0
!
!    Input arguments:
!  n             : the neutrino flavor index
!  j             : radial zone index
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  agrjmh(j)     : lapse function at zone j-1/2 before x-hydro step
!  agrajmh(j)    : lapse function at zone j-1/2 after x-hydro step
!  unub(j,k)     : zone-edged energy for mass zone j, energy zone k
!  nnugp(n)      : number of neutrino energy groups for n-neutrinos
!
!    Output arguments (common):
!  v_e(j,k)      : neutrino advection velocity through the energy grid
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half

USE e_advct_module, ONLY : adot, adot_a, v_e, v_e1, agrjmh, agrajmh
USE nu_energy_grid_module, ONLY : nnugp, unubi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: j             ! radial zone index
INTEGER, INTENT(in)               :: n             ! neutrino flavor index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                           :: k             ! neutrino energy index

REAL(KIND=double)                 :: unub          ! zone-edged energy
REAL(KIND=double)                 :: agrjmh_1      ! 1/agrjmh

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute adot(j), the time differences of agrjmh
!  Compute adot_a(j) the time differences  of agrjmh divided by its 
!   value.
!-----------------------------------------------------------------------

adot(j)           = agrajmh(j) - agrjmh(j)
adot_a(j)         = adot(j)/( half * ( agrajmh(j) + agrjmh(j) ) )

!-----------------------------------------------------------------------
!  Compute the energy zone-edged value of the advection displacement
!   v_e1.
!-----------------------------------------------------------------------

v_e1(j)           = adot_a(j)

!-----------------------------------------------------------------------
!  Compute the energy zone-edged value of the advection displacement.
!-----------------------------------------------------------------------

agrjmh_1          = 1.d0/agrjmh(j)
DO k = 2,nnugp(n)
  unub            = unubi(k) * agrjmh_1
  v_e(j,k)        = half * ( agrajmh(j) + agrjmh(j) ) * unub * v_e1(j)
END DO

v_e(j,1)          = zero
v_e(j,nnugp(n)+1) = zero

RETURN
END SUBROUTINE e_advct_agr_x_v

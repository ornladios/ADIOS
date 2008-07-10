SUBROUTINE e_advct_x_v( j, n )
!-----------------------------------------------------------------------
!
!    File:         e_advct_x_v
!    Module:       e_advct_x_v
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/26/01
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
!  r(j)          : radius of zone j before x-hydro step (cm)
!  ra(j)         : radius of zone j after x-hydro step (cm)
!  rho(j)        : matter density of zone j before x-hydro step (g/cm**3)
!  rhoa(j)       : matter density of zone j after x-hydro step (g/cm**3)
!  agrjmh(j)     : lapse function at zone j-1/2 before x-hydro step
!  agrajmh(j)    : lapse function at zone j-1/2 after x-hydro step
!  unub(j,k)     : zone-edged energy for mass zone j, energy zone k
!  nnugp(n)      : number of neutrino energy groups for n-neutrinos
!  Ef(j,k,n)     : Eddington factor as a function of j,k,n
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
USE numerical_module, ONLY : zero, half, third, epsilon

USE e_advct_module, ONLY : adot, ddot, rdot, adot_a, ddot_d, rdot_r, &
& bdot_b, v_e, v_e1, v_e2, rjmh, rajmh, r, ra, rho, rhoa, agrjmh, agrajmh, Ef
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
!  Interpolate r & ra to the zone centers
!-----------------------------------------------------------------------

rjmh(j)           = half * ( ( r(j)  + epsilon )**3 + ( r(j-1)  + epsilon )**3 )**third
rajmh(j)          = half * ( ( ra(j) + epsilon )**3 + ( ra(j-1) + epsilon )**3 )**third

!-----------------------------------------------------------------------
!  Compute  adot(j), ddot(j), rdot(j), the time differences of agrjmh,
!   rho, and rjmh.
!  Compute  adot_a(j), bdot_b(j), ddot_d(j), and rdot_r(j), the time
!   differences  of agrjmh, b, rho, and rjmh divided by their respective
!   values.
!-----------------------------------------------------------------------

adot(j)           = agrajmh(j) - agrjmh(j)
ddot(j)           = rhoa(j)    - rho(j)
rdot(j)           = rajmh(j)   - rjmh(j)

adot_a(j)         = adot(j)/( half * ( agrajmh(j) + agrjmh(j) ) )
ddot_d(j)         = ddot(j)/( half * ( rhoa(j)    + rho(j)    ) )
rdot_r(j)         = rdot(j)/( half * ( rajmh(j)   + rjmh(j)   ) )
bdot_b(j)         = - 2.d0 * rdot_r(j) - ddot_d(j)

!-----------------------------------------------------------------------
!  Compute the energy zone-edged value of the advection displacement
!   v_e1 and v_e2.
!-----------------------------------------------------------------------

v_e1(j)           = adot_a(j) - rdot_r(j)
v_e2(j)           = bdot_b(j) - rdot_r(j)

!-----------------------------------------------------------------------
!  Compute the energy zone-edged value of the advection displacement.
!-----------------------------------------------------------------------

agrjmh_1          = 1.d0/agrjmh(j)
DO k = 2,nnugp(n)
  unub            = unubi(k) * agrjmh_1
  v_e(j,k)        = half * ( agrajmh(j) + agrjmh(j) ) * unub * ( v_e1(j) - Ef(j,k,n) * v_e2(j) )
END DO

v_e(j,1)          = zero
v_e(j,nnugp(n)+1) = zero

RETURN
END SUBROUTINE e_advct_x_v

SUBROUTINE e_advct_y_v( j, n, i_radial )
!-----------------------------------------------------------------------
!
!    File:         e_advct_y_v
!    Module:       e_advct_y_v
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/06/05
!
!    Purpose:
!      To compute the neutrino advection velocity through the
!       energy grid.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  n             : the neutrino flavor index
!  j             : y (angular) zone index
!  i_radial      : the unshifted radial zone (angular ray) corresponding to ji_ray, jk_ray
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  y(j)          : y-coordinate of zone j before y-hydro step (cm)
!  ya(j)         : y-coordinate of zone j after y-hydro step (cm)
!  rho(j)        : matter density of zone j before y-hydro step (g/cm**3)
!  rhoa(j)       : matter density of zone j after y-hydro step (g/cm**3)
!  unub(j,k)     : zone-edged energy for mass zone j, energy zone k (MeV)
!  nnugp(n)      : number of neutrino energy groups for n-neutrinos
!
!    Output arguments (common):
!  v_e(j,k)      : neutrino advection velocity through the energy grid
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, nu_energy_grid_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half, third, epsilon

USE e_advct_module, ONLY : adot, ddot, rdot, adot_a, ddot_d, rdot_r, &
& bdot_b, v_e, v_e1, v_e2, rjmh, yjmh, yajmh, y, ya, rho, rhoa
USE nu_energy_grid_module, ONLY : nnugp, unubi
USE parallel_module, ONLY : myid

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: j             ! y (angular) zone index
INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: i_radial      ! the unshifted radial zone corresponding to ji_ray, jk_ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                           :: k             ! neutrino energy index

REAL(KIND=double)                 :: Ef            ! Eddington factor
REAL(KIND=double)                 :: unub          ! zone-edged energy

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set Eddington factor to 1/3 to approximate small mean free paths
!-----------------------------------------------------------------------

Ef                = third

!-----------------------------------------------------------------------
!  Physical values of y and ya.
!-----------------------------------------------------------------------

yjmh(j)           = ( y(j)  + y(j+1)  ) * rjmh(i_radial)
yajmh(j)          = ( ya(j) + ya(j+1) ) * rjmh(i_radial)

!-----------------------------------------------------------------------
!  Compute  adot(j), ddot(j), rdot(j), the time differences of agrjmh,
!   rho, and yjmh.
!  Compute  adot_a(j), bdot_b(j), ddot_d(j), and rdot_r(j), the time
!   differences  of agrjmh, b, rho, and yjmh divided by their respective
!   values.
!-----------------------------------------------------------------------

ddot(j)           = rhoa(j)    - rho(j)
rdot(j)           = yajmh(j)   - yjmh(j)

ddot_d(j)         = ddot(j)/( half * ( rhoa(j)    + rho(j)    ) )
rdot_r(j)         = rdot(j)/( half * ( yajmh(j)   + yjmh(j)   ) )
bdot_b(j)         = - 2.d0 * rdot_r(j) - ddot_d(j)

!-----------------------------------------------------------------------
!  Compute the energy zone-edged value of the advection displacement
!   v_e1 and v_e2.
!-----------------------------------------------------------------------

v_e1(j)           = - rdot_r(j)
v_e2(j)           = bdot_b(j) - rdot_r(j)

!-----------------------------------------------------------------------
!  Compute the energy zone-edged value of the advection displacement.
!-----------------------------------------------------------------------

DO k = 2,nnugp(n)
  unub            = unubi(k)
  v_e(j,k)        = unub * ( v_e1(j) - Ef * v_e2(j) )
END DO

v_e(j,1)          = zero
v_e(j,nnugp(n)+1) = zero

RETURN
END SUBROUTINE e_advct_y_v

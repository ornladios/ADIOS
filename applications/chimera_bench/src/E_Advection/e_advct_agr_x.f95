SUBROUTINE e_advct_agr_x( jr_min, jr_max, nn )
!-----------------------------------------------------------------------
!
!    File:         e_advct_agr_x
!    Module:       e_advct_agr_x
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/19/06
!
!    Purpose:
!      To compute the neutrino advection velocity through the
!       energy grid as a result of the x-hydro step.
!
!    Subprograms called:
!  e_advct_agr_x_v  : calculates neutrino "energy velocities" at the energy zone edges
!  e_advct_bc       : computes psi0 boundary values
!  e_advct_vol      : calculates the energy space volumes
!  e_advct_evol     : updates the energy zone boundaries and psi0 (Lagrangian step)
!  e_advct_remap    : remaps the energy zone boundaries back to the Eulerian grid
!
!    Input arguments:
!  jr_min           : inner zone for which nu energy advection is to be calculated
!  jr_max           : outer zone for which nu energy advection is to be calculated
!  nn               : the neutrino flavor index
!
!    Output arguments:
!        none
!
!    Input arguments (common)
!  rho(j)           : matter density of zone j before the x-hydro step (g/cm**3)
!  rhoa(j)          : matter density of zone j after the x-hydro step (g/cm**3)
!  unubi            : zone-edged energy for mass zone j, energy zone k at infity
!  dunui            : unub(j,k+1) - unub(j,k) at infity
!  nnugp(n)         : number of neutrino energy groups for n-neutrinos
!
!    Output arguments (common):
!  nmin             : number of first real zone (=1+6)
!  nmax             : number of last  real zone (=nnugp(n)+6)
!  xk(k)            : energy grid zone edge locations
!  dk(k)            : xk(k+1) - xk(k)
!  xk_0(k)          : energy grid zone edge locations prior to lagrangian update
!  dk_0(k)          : dk(k) prior to lagrangian update
!  xk_1(k)          : energy grid zone edge locations prior at m+1
!  dk_1(k)          : dk(k) at m+1
!  psi(k)           : neutrino occupation number
!  v_e(j,k)         : neutrino advection velocity through the energy grid
!  psi0_a           : zero angular moments of the neutrino occupation number
!
!    Include files:
!  numerical_module
!  e_advct_module, evh1_global, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero
	
USE e_advct_module, ONLY : nmin, nmax, xk_0, dk_0, xk_1, dk_1, xk, dk, vk, &
& psi, v_e, psi0, psi0_a, agrjmh, agrajmh
USE evh1_global, ONLY : ngeomx
USE nu_energy_grid_module, ONLY : nnugp, unui, dunui, unubi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nn            ! neutrino flavor index
INTEGER, INTENT(in)               :: jr_min        ! minimum radial zone
INTEGER, INTENT(in)               :: jr_max        ! maximum radial zone

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                           :: j             ! radial zone index
INTEGER                           :: ntot          ! nnugp(nn) + 12

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Add six ghost zones to 1D sweep array
!-----------------------------------------------------------------------

nmin                     = 7
nmax                     = nnugp(nn) + 6
ntot                     = nnugp(nn) + 12
 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                  ||||| Begin radial loop |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
 
!-----------------------------------------------------------------------
!
!            \\\\\ COMPUTE NEUTRINO ENERGY ADVECTION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

  psi (1:ntot)           = zero
  xk_0(1:ntot)           = zero
  dk_0(1:ntot)           = zero
  xk_1(1:ntot)           = zero
  dk_1(1:ntot)           = zero
  xk  (1:ntot)           = zero
  dk  (1:ntot)           = zero
  vk  (1:ntot)           = zero

!-----------------------------------------------------------------------
!  Compute the energy zone edge velocities
!-----------------------------------------------------------------------

  CALL e_advct_agr_x_v(j,nn)
 
!-----------------------------------------------------------------------
!  Put variables into 1D arrays, padding with 6 ghost zones
!-----------------------------------------------------------------------

  psi (7:nnugp(nn)+6)    = psi0(j,1:nnugp(nn),nn)/agrjmh(j)**3
  xk_0(7:nnugp(nn)+6)    = unubi(1:nnugp(nn))
  dk_0(7:nnugp(nn)+6)    = dunui(1:nnugp(nn))
  xk_1(7:nnugp(nn)+6)    = unubi(1:nnugp(nn))
  dk_1(7:nnugp(nn)+6)    = dunui(1:nnugp(nn))
  xk  (7:nnugp(nn)+6)    = unubi(1:nnugp(nn))
  dk  (7:nnugp(nn)+6)    = dunui(1:nnugp(nn))
  vk  (7:nnugp(nn)+6)    = v_e (j,1:nnugp(nn))

  vk (nnugp(nn)+7)       = v_e (j,nnugp(nn)+1)
 
!-----------------------------------------------------------------------
!  Impose the boundary conditions
!-----------------------------------------------------------------------

  CALL e_advct_bc
 
!-----------------------------------------------------------------------
!  Calculate the energy space volume
!-----------------------------------------------------------------------

  CALL e_advct_vol
 
!-----------------------------------------------------------------------
!  Update the energy zone boundaries and psi0
!-----------------------------------------------------------------------

  CALL e_advct_evol
 
!-----------------------------------------------------------------------
!  Remap the energy zone boundaries back to the Eulerian grid
!-----------------------------------------------------------------------

  CALL e_advct_remap( ngeomx )
 
!-----------------------------------------------------------------------
!  Store updated psi0 in psi0_a
!-----------------------------------------------------------------------

  psi0_a(j,1:nnugp(nn),nn) = psi(7:nnugp(nn)+6) * agrajmh(j)**3

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                   ||||| End radial loop |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END DO

RETURN
END SUBROUTINE e_advct_agr_x

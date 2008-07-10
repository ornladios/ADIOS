SUBROUTINE e_advct_z( kmin, kmax, nn, ki_ray, kj_ray, i_radial )
!-----------------------------------------------------------------------
!
!    File:         e_advct_z
!    Module:       e_advct_z
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/07
!
!    Purpose:
!      To compute the neutrino advection velocity through the
!       energy grid as a result of the z-hydro step.
!
!    Subprograms called:
!  e_advct_z_v      : calculates neutrino "energy velocities" at the energy zone edges
!  e_advct_bc       : computes psi0 boundary values
!  e_advct_vol      : calculates the energy space volumes
!  e_advct_evol     : updates the energy zone boundaries and psi0 (Lagrangian step)
!  e_advct_remap    : remaps the energy zone boundaries back to the Eulerian grid
!
!    Input arguments:
!  kmin             : minimum z-array index
!  kmax             : maximum z-array index
!  nn               : the neutrino flavor index
!  ki_ray           : x (radial) index of a specific z (azimuthal) ray
!  kj_ray           : y (angular) index of a specific z (azimuthal) ray
!  i_radial         : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!
!    Output arguments:
!        none
!
!    Input arguments (common)
!  rho(j)           : matter density of zone j before the z-hydro step (g/cm**3)
!  rhoa(j)          : matter density of zone j after the z-hydro step (g/cm**3)
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
!  unujv(j,n,i_ray) : energy transferred to n-type neutrinos by energy advection (ergs)
!  unucrv(nj_ray)   : total energy transferred to n-type neutrinos by energy advection (ergs)
!  dunujvdt(j,n)    : rate of energy transferred to n-type neutrinos by energy advection
!  dndt_v(j,k,n)    : net rate of n-neutrino production by energy advection
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  e_advct_module, evh1_global, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half, one, third, epsilon, ncoef, ecoef
USE physcnst_module, ONLY : cvel
	
USE e_advct_module, ONLY : nmin, nmax, xk_0, dk_0, xk_1, dk_1, xk, dk, vk, &
& psi, v_e, u, ivc_y, dtnph, psi0, psi0_a, rho, rhoa, dmrst, ncoefa, ncoefaa, &
& ecoefa, ecoefaa, unujv, unucrv, dunujvdt, dndt_v, psi1
USE evh1_global, ONLY : ngeomx
USE nu_energy_grid_module, ONLY : nnugp, unui, dunui, unubi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nn            ! neutrino flavor index
INTEGER, INTENT(in)               :: kmin          ! minimum radial zone
INTEGER, INTENT(in)               :: kmax          ! maximum radial zone
INTEGER, INTENT(in)               :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)               :: kj_ray        ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)               :: i_radial      ! the unshifted radial zone corresponding to ki_ray, kj_ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                           :: j             ! radial zone index
INTEGER                           :: ntot          ! nnugp(nn) + 12
INTEGER                           :: it            ! iteration index

REAL(KIND=double)                 :: d_rho         ! ( rhoa/rho )**1/3
REAL(KIND=double)                 :: u_nu_i        ! initial energy of n-neutrinos
REAL(KIND=double)                 :: u_nu_ideal    ! final energy of n-neutrinos given by rho**1/3 law
REAL(KIND=double)                 :: u_nu_f        ! final energy of n-neutrinos
REAL(KIND=double)                 :: unuvt         ! change in neutrino energy/gram due to neutrino energy advection
REAL(KIND=double)                 :: v_scale       ! ratio of rho**1/3 law to actual energy advance

REAL(KIND=double), PARAMETER      :: rho_high = 3.d9

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Add six ghost zones to 1D sweep array
!-----------------------------------------------------------------------

nmin                      = 7
nmax                      = nnugp(nn) + 6
ntot                      = nnugp(nn) + 12
 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                ||||| Begin z (azimuthal) loop |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO j = kmin,kmax
  
!-----------------------------------------------------------------------
!
!            \\\\\ COMPUTE NEUTRINO ENERGY ADVECTION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

  unuvt                   = zero

  psi (1:ntot)            = zero
  xk_0(1:ntot)            = zero
  dk_0(1:ntot)            = zero
  xk_1(1:ntot)            = zero
  dk_1(1:ntot)            = zero
  xk  (1:ntot)            = zero
  dk  (1:ntot)            = zero
  vk  (1:ntot)            = zero

  ncoefa (j,1:nnugp(nn))  = ncoef * unui(1:nnugp(nn))**2 * dunui(1:nnugp(nn))
  ecoefa (j,1:nnugp(nn))  = ecoef * unui(1:nnugp(nn))**3 * dunui(1:nnugp(nn))
  ncoefaa(j,1:nnugp(nn))  = ncoef * unui(1:nnugp(nn))**2 * dunui(1:nnugp(nn))
  ecoefaa(j,1:nnugp(nn))  = ecoef * unui(1:nnugp(nn))**3 * dunui(1:nnugp(nn))

!-----------------------------------------------------------------------
!  Compute the energy zone edge velocities
!-----------------------------------------------------------------------

  CALL e_advct_z_v( j, nn, i_radial )
 
!-----------------------------------------------------------------------
!  Put variables into 1D arrays, padding with 6 ghost zones
!-----------------------------------------------------------------------

  psi (7:nnugp(nn)+6)     = psi0(j,1:nnugp(nn),nn)/rho(j)
  xk_0(7:nnugp(nn)+6)     = unubi(1:nnugp(nn))
  dk_0(7:nnugp(nn)+6)     = dunui(1:nnugp(nn))
  xk_1(7:nnugp(nn)+6)     = unubi(1:nnugp(nn))
  dk_1(7:nnugp(nn)+6)     = dunui(1:nnugp(nn))
  xk  (7:nnugp(nn)+6)     = unubi(1:nnugp(nn))
  dk  (7:nnugp(nn)+6)     = dunui(1:nnugp(nn))
  vk  (7:nnugp(nn)+6)     = v_e (j,1:nnugp(nn))

  vk (nnugp(nn)+7)        = v_e (j,nnugp(nn)+1)
 
!-----------------------------------------------------------------------
!  Impose the boundary conditions
!-----------------------------------------------------------------------

  CALL e_advct_bc
 
!-----------------------------------------------------------------------
!  Update the energy zone boundaries and psi0
!-----------------------------------------------------------------------

  CALL e_advct_evol
 
!-----------------------------------------------------------------------
!  Remap the energy zone boundaries back to the Eulerian grid
!-----------------------------------------------------------------------

  CALL e_advct_remap( ngeomx )
  
!-----------------------------------------------------------------------
!
!         \\\\\           IF RHO > RHO_HIGH,             /////
!         \\\\\ NU_ENERGY DENSITY MUST SCALE AS RHO**4/3 /////
!
!-----------------------------------------------------------------------

  IF ( rho(j) >= rho_high  .and.  rhoa(j) /= rho(j) ) THEN

!-----------------------------------------------------------------------
!        ||||| Iterate twice to get achieve proper scaling |||||
!-----------------------------------------------------------------------
    
    DO it = 1, 2

!-----------------------------------------------------------------------
!  Compute density increment
!-----------------------------------------------------------------------

      d_rho               = ( rhoa(j)/rho(j) )**third
 
!-----------------------------------------------------------------------
!  Compute initial and final energy of n-neutrinos, and the energy
!   increment
!-----------------------------------------------------------------------

      u_nu_i              = zero
      u_nu_f              = zero

      u_nu_i              = SUM( ecoefa(j,1:nnugp(nn)) * psi0(j,1:nnugp(nn),nn) )
      u_nu_f              = SUM( ecoefaa(j,1:nnugp(nn)) * psi(7:nnugp(nn)+6) ) * rhoa(j)

      u_nu_i              = u_nu_i/rho(j)
      u_nu_f              = u_nu_f/rhoa(j)
      u_nu_ideal          = d_rho * u_nu_i
 
!-----------------------------------------------------------------------
!  Rescale the energy zone edge velocities and reinitialize variables
!-----------------------------------------------------------------------

      psi (1:ntot)        = zero
      xk_0(1:ntot)        = zero
      dk_0(1:ntot)        = zero
      xk  (1:ntot)        = zero
      dk  (1:ntot)        = zero

      v_scale             = ( u_nu_ideal - u_nu_i )/( u_nu_f - u_nu_i + epsilon )
      v_scale             = DMAX1( DMIN1( v_scale, 1.0d+01 ), 1.0d-01 )

      psi (7:nnugp(nn)+6) = psi0(j,1:nnugp(nn),nn)/rho(j)
      xk_0(7:nnugp(nn)+6) = unubi(1:nnugp(nn))
      dk_0(7:nnugp(nn)+6) = dunui(1:nnugp(nn))
      xk  (7:nnugp(nn)+6) = unubi(1:nnugp(nn))
      dk  (7:nnugp(nn)+6) = dunui(1:nnugp(nn))
      vk  (7:nnugp(nn)+6) = v_scale * vk(7:nnugp(nn)+6)
 
!-----------------------------------------------------------------------
!  Impose the boundary conditions
!-----------------------------------------------------------------------

      CALL e_advct_bc
 
!-----------------------------------------------------------------------
!  Update the energy zone boundaries and psi0
!-----------------------------------------------------------------------

      CALL e_advct_evol
 
!-----------------------------------------------------------------------
!  Remap the energy zone boundaries back to the Eulerian grid
!-----------------------------------------------------------------------

      CALL e_advct_remap( ngeomx )

!-----------------------------------------------------------------------
!                     ||||| End iteration |||||
!-----------------------------------------------------------------------

    END DO ! it = 1, 2

  END IF ! rho(j) >= rho_high  .and.  rhoa(j) /= rho(j)
 
!-----------------------------------------------------------------------
!  Store updated psi0 in psi0_a
!-----------------------------------------------------------------------

  psi0_a(j,1:nnugp(nn),nn) = psi(7:nnugp(nn)+6) * rhoa(j)

!-----------------------------------------------------------------------
!  dndt_v(j,k,n) : net rate of n-neutrino production by energy advection
!                  in zones j, k and neutrino type n 
!                  (cm^{-3} s^{-1} MeV^{-1})
!-----------------------------------------------------------------------

!  DO k = 1,nnugp(nn)
!    dndt_v(j,k,nn,kj_ray,ki_ray) = ncoefaa(j,k) * psi0_a(j,k,nn) - ncoefa(j,k) * psi0(j,k,nn)
!  END DO

!-----------------------------------------------------------------------
!  unujv(j,n)      : energy transferred to n-type neutrinos by advection
!                     in zone j (ergs)
!  unucrv(n,kj_ray,ki_ray) : total energy transferred to n-type neutrinos by
!                     advection (ergs)
!  dunujvdt(j,n)   : rate of energy transferred to n-type neutrinos by
!                     advection in zone j (ergs s^{-1} g^{-1})
!-----------------------------------------------------------------------

!  DO k = 1,nnugp(nn)
!    flx_fact         = DMAX1( DMIN1( half * ( psi1(j,k,nn) + psi1(j-1,k,nn) )/( psi0(j,k,nn) + epsilon ), one ), zero )
!    unuvtt           = ecoefaa(j,k) * psi0_a(j,k,nn) * ( one + flx_fact * u(j)/cvel )/rhoa(j) & 
!&                    - ecoefa(j,k)  * psi0(j,k,nn)   * ( one + flx_fact * u(j)/cvel )/rho(j)
!    unuvt            = unuvt + unuvtt
!  END DO

!  unujv(j,nn,kj_ray,ki_ray)  = unujv(j,nn,kj_ray,ki_ray) + unuvt * dmrst(j)
!  unucrv(nn,kj_ray,ki_ray)   = unucrv(nn,kj_ray,ki_ray) + unuvt * dmrst(j)
!  dunujvdt(j,nn,kj_ray,ki_ray) = unuvt/dtnph

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                 ||||| End z (azimuthal) loop |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END DO

RETURN
END SUBROUTINE e_advct_z

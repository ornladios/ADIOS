SUBROUTINE psi1_cal( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& rho, t, ye, r, rstmss, u, psi0p, psi1p, nx, nez, nnu, it )
!-----------------------------------------------------------------------
!
!    File:         psi1_cal
!    Module:       psi1_cal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/12/00
!
!    Purpose:
!      To compute psi1 from psi0 and the diffusion coefficients.
!
!    Subprograms called:
!  abemrate     : computes emission and absorption rates for mfp's
!  bremrate     : computes contribution to RHS from nucleon-nucleon bremsstrahlung
!  pairrate     : computes contribution to RHS from electron-positron pair annihilation
!  scterate     : computes contribution to RHS from neutrino-electron scattering
!  sctirate     : computes isoenergetic scattering rates for mfp's
!  sctnnrate    : computes contribution to RHS from neutrino-nucleon inelastic scattering
!  sctnArate    : computes contribution to RHS from neutrino-nucleous inelastic scattering
!  sctnrate     : computes contribution to RHS from neutrino-nucleon elastic scattering
!  mfp_cal      : computes inverse mfp's
!  pre_trans    : computes quantities needed to compute diffusion coefficients
!  nu_sphere    : computes neutrinospheres, needed to compute diffusion coefficients
!  diffcr       : computes diffusion coefficients if first if Case 1
!  diff         : computes diffusion coefficients for the current time step
!  gamgr_nu_cal : computes relativistic gammas for computing proper lengths
!
!    Input arguments:
!
!  jr_min       : inner x-array index
!  jr_max       : outer x-array index
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  n            : neutrino type
!  k            : neutrino energy group index
!  rho          : density (g cm^{-3})
!  t            : temperature (K)
!  ye           : electron fraction
!  r            : radius (cm)
!  rstmss       : enclosed rest mass (g)
!  u            : radial velocity (cm s^{-1}
!  psi0p        : zeroth angular moment of theneutrino distribution function
!  nx           : x-array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  it           : iteration count
!
!    Output arguments:
!  psi1p        : first angular moment of theneutrino distribution function
!
!    Input arguments (common):
!  dc(j,k,n)    : diffusion coefficient
!
!    Include files:
!  kind_module, numerical_module
!  nu_dist_module, nu_energy_grid_module, prb_cntl_module,
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half

USE nu_dist_module, ONLY : dc, dcr, psi0, psi1, drjmh_inv, j_sphere, &
& r_sphere, d_sphere, t_sphere, m_sphere
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY : idiff 

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nx            ! array dimension
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: it            ! iteration count

REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: rho    ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: t      ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: ye     ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: r      ! radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: rstmss ! enclosed rest mass (g)
REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: u      ! x-velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu)    :: psi0p  ! zero moment of the NDF

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu) :: psi1p  ! first moment of the NDF

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index

REAL(KIND=double)                :: psi_ratio     ! ratio of psi0 in outer ghost to outer psi0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer neutrino distribution functions
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max
      psi0(j,k,n)  = psi0p(j,k,n)
      psi1(j,k,n)  = psi1p(j,k,n)
    END DO ! j = jr_min,jr_max
  END DO ! k = 1,nnugp(n)
END DO ! n = 1,nnu

!-----------------------------------------------------------------------
!  Absorption and emission inverse mean free paths
!-----------------------------------------------------------------------

CALL abemrate( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )

!-----------------------------------------------------------------------

DO n = 1,nnu

  IF ( nnugp(n) /= 0 ) THEN

!-----------------------------------------------------------------------
!  Isoenergetic scattering inverse mean free paths
!-----------------------------------------------------------------------

    CALL sctirate( jr_min, jr_max, ij_ray, ik_ray, n, rho, t, ye, nx )

!-----------------------------------------------------------------------
!  e-neutrino-electron scattering
!-----------------------------------------------------------------------

    CALL scterate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )

!-----------------------------------------------------------------------
!  e-neutrino production from electron-positron pair annihilation
!-----------------------------------------------------------------------

    CALL pairrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, r, nx )

!-----------------------------------------------------------------------
!  e-neutrino production from nucleon-nucleon pair bremsstrahlung
!   annihilation
!-----------------------------------------------------------------------

    CALL bremrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, r, nx )

!-----------------------------------------------------------------------
!  e-neutrino-nucleon elastic scatteringn
!-----------------------------------------------------------------------

    CALL sctnrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )

!-----------------------------------------------------------------------
!  e-neutrino-nucleon inelastic scattering
!-----------------------------------------------------------------------

    CALL sctnnrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )

!-----------------------------------------------------------------------
!  e-eutrino-nucleus inelastic scattering
!-----------------------------------------------------------------------

    CALL sctnArate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )

  END IF ! nnugp(n) /= 0

END DO ! n = 1,nnu

!-----------------------------------------------------------------------
!  Inverse mean free paths
!-----------------------------------------------------------------------

CALL mfp_cal( jr_min, jr_max, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Neutrinospheres
!-----------------------------------------------------------------------

CALL pre_trans( jr_min, jr_max, rho, r, nx, nnu )
CALL nu_sphere( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& r, rho, t, rstmss, nx, nez, nnu, j_sphere, r_sphere, d_sphere, t_sphere, &
& m_sphere )

!-----------------------------------------------------------------------
!  Call diffcr if first = true
!-----------------------------------------------------------------------

IF ( first ) THEN
  CALL diffcr( jr_min, jr_max, ij_ray, ik_ray, r, u, nx )
  first            = .false.
END IF

!-----------------------------------------------------------------------
!  Outer boundary condition
!-----------------------------------------------------------------------

CALL psi_bd( jr_max, 1, 1, r, nx, psi_ratio )

!-----------------------------------------------------------------------
!
!                  \\\\\ DIFFUSION COEFFICIENT ////
!
!        The diffusion coefficients dc(j,k,n) are updated by
!         subroutine diffc.
!
!    idiff = 0: dc(j,k,n) computed at each iteration
!    idiff = 1: dc(j,k,n) is the average of the current and preceding timestep value
!    idiff = 2: dc(j,k,n) = 0 for j = jr_max, dc(j,k,n) as in idiff = 0 for j /= jr_max
!    idiff = 3: dc(j,k,n) = 0 for j = jr_max, dc(j,k,n) as in idiff = 1 for j /= jr_max
!    idiff = 4: dc(j,k,n) = 0 for all j
!    idiff = 5: dc(j,k,n) is computed only on the first iteration
!    idiff = 6: dc(j,k,n) = 0 for j = jr_max, computed only on the first iteration for j /= jr_max
!
!    dcr(j,k,n) is obtained from the preceding dc(j,k,n) in subroutine nu_trans
!-----------------------------------------------------------------------

SELECT CASE (idiff)

  CASE(0)

    CALL diffc( jr_min, jr_max, ij_ray, ik_ray, r, u, nx )

  CASE(1)

    CALL diffc( jr_min, jr_max, ij_ray, ik_ray, r, u, nx )
    DO n = 1,nnu
      DO k = 1,nnugp(n)
        DO j = jr_min,jr_max
          dc(j,k,n)  = half * ( dc(j,k,n) + dcr(j,k,n) )
        END DO ! j = jr_min,jr_max
      END DO ! k = 1,nnugp(n)
    END DO ! n = 1,nnu

  CASE(2)

    CALL diffc( jr_min, jr_max-1, ij_ray, ik_ray, r, u, nx )
    DO n = 1,nnu
      DO k = 1,nnugp(n)
        dc(jr_max,k,n) = zero
      END DO ! k = 1,nnugp(n)
    END DO ! n = 1,nnu

  CASE(3)

    CALL diffc( jr_min, jr_max-1, ij_ray, ik_ray, r, u, nx )
    DO n = 1,nnu
      DO k = 1,nnugp(n)
        DO j = jr_min,jr_max-1
          dc(j,k,n)  = half * ( dc(j,k,n) + dcr(j,k,n) )
        END DO ! j = jr_min,jr_max-1
      END DO ! k = 1,nnugp(n)
    END DO ! n = 1,nnu

    DO n = 1,nnu
      DO k = 1,nnugp(n)
        dc(jr_max,k,n) = zero
      END DO ! k = 1,nnugp(n)
    END DO ! n = 1,nnu

  CASE(4)

    DO n = 1,nnu
      DO k = 1,nnugp(n)
        DO j = jr_min,jr_max
          dc(j,k,n)  = zero
        END DO ! j = jr_min,jr_max
      END DO ! k = 1,nnugp(n)
    END DO ! n = 1,nnu

  CASE(5)

    IF ( it == 1 ) THEN
      CALL diffc( jr_min, jr_max, ij_ray, ik_ray, r, u, nx )
    END IF ! it == 1

  CASE(6)

    IF ( it == 1 ) THEN

      CALL diffc( jr_min, jr_max, ij_ray, ik_ray, r, u, nx )

      DO n = 1,nnu
        DO k = 1,nnugp(n)
          dc(jr_max,k,n) = zero
        END DO ! k = 1,nnugp(n)
      END DO ! n = 1,nnu

    END IF ! it == 1

END SELECT

!-----------------------------------------------------------------------
!
!                           \\\\\ PSI1 ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Get ghost values of psi0
!-----------------------------------------------------------------------

  IF ( nnugp(1) /= 0 ) THEN
    DO k = 1,nnugp(1)
      CALL psi_bd( jr_max, k, 1, r, nx, psi_ratio )
      psi0(jr_max+1,k,1) = psi_ratio * psi0(jr_max,k,1)
    END DO ! k = 1,nnugp(1)
  END IF ! nnugp(1) /= 0
    
  IF ( nnugp(2) /= 0 ) THEN
    DO k = 1,nnugp(2)
      CALL psi_bd( jr_max, k, 2, r, nx, psi_ratio )
      psi0(jr_max+1,k,2) = psi_ratio * psi0(jr_max,k,2)
    END DO ! k = 1,nnugp(2)
  END IF ! nnugp(2) /= 0
    
  IF ( nnugp(3) /= 0 ) THEN
    DO k = 1,nnugp(3)
      CALL psi_bd( jr_max, k, 3, r, nx, psi_ratio )
      psi0(jr_max+1,k,3) = psi_ratio * psi0(jr_max,k,3)
    END DO ! k = 1,nnugp(3)
  END IF ! nnugp(3) /= 0
    
  IF ( nnugp(4) /= 0 ) THEN
    DO k = 1,nnugp(4)
      CALL psi_bd( jr_max, k, 4, r, nx, psi_ratio )
      psi0(jr_max+1,k,4) = psi_ratio * psi0(jr_max,k,4)
    END DO ! k = 1,nnugp(4)
  END IF ! nnugp(4) /= 0

!-----------------------------------------------------------------------
!  Now compute psi1
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max
      psi1(j,k,n)    = - dc(j,k,n) * ( psi0(j+1,k,n) - psi0(j,k,n) ) * drjmh_inv(j)
    END DO ! j = jr_min,jr_max
  END DO ! k = 1,nnugp(n)
END DO ! n = 1,nnu

IF ( idiff == 4  .or.  idiff == 6 ) THEN
  DO n = 1,nnu
    DO k = 1,nnugp(n)
      psi1(jr_max,k,n) = zero
    END DO ! j = jr_min,jr_max
  END DO ! n = 1,nnu
END IF ! n = 1,nnu

!-----------------------------------------------------------------------
!  Transfer neutrino distribution functions out
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max
      psi1p(j,k,n) = psi1(j,k,n)
    END DO ! j = jr_min,jr_max
  END DO ! k = 1,nnugp(n)
END DO ! n = 1,nnu

RETURN
END SUBROUTINE psi1_cal

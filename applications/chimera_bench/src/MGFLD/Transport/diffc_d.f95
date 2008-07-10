SUBROUTINE diffc_d( jr_min, jr_max, ij_ray, ik_ray, r, psi0, p_dpsi0, &
& dc_jph, dc_jmh, nx, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         diffc_d
!    Module:       diffc_d
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/19/00
!
!    Purpose:
!      To compute the diffusion coefficient, dc
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min         : inner radial zone for computing the diffusion coefficient.
!  jr_max         : outer radial zone for computing the diffusion coefficient.
!  ij_ray         : j-index of a radial ray
!  ik_ray         : k-index of a radial ray
!  r              : radius (cm)
!  psi0(j,k,n)    : neutrino distribution function
!  p_dpsi0(j,k,n) : perturbd neutrino distribution function
!  nx             : x-array extent
!  nez            : neutrino energy array extent
!  nnu            : neutrino flavor array extent
!
!    Output arguments:
!  dc_jph(j,k,n) : diffusion coefficient when psi0(j+1,k,n) is perturbed
!  dc_jmh(j,k,n) : diffusion coefficient when psi0(j  ,k,n) is perturbed
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half, one, third, epsilon
USE physcnst_module, ONLY : cvel, g

USE cycle_module
USE edit_module, ONLY : nlog
USE mdl_cnfg_module, ONLY : u, grvmss
USE nu_dist_module, ONLY : tmfp, r_sphere, m_sphere, rjmh, drjmh_inv, &
& dc, fluxlm, tmfp_j, j_sphere
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY : ireltrns

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: nx            ! radial array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)            :: r        ! radius (cm)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu)   :: dc_jph   ! diffusion coefficient when psi0(j+1,k,n) is perturbed
REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu)   :: dc_jmh   ! diffusion coefficient when psi0(j  ,k,n) is perturbed

!-----------------------------------------------------------------------
!        Input - Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu) :: psi0     ! neutrino distribution function
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu) :: p_dpsi0  ! perturbd neutrino distribution function

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status

REAL(KIND=double)                :: psi0j         ! psi0 interpolated to zone edge
REAL(KIND=double)                :: dpsidr        ! gradient in psi0
REAL(KIND=double)                :: guiuv         ! grad psi0/psi0
REAL(KIND=double)                :: tmfp_jj       ! mean free path interpolater to zone edge
REAL(KIND=double)                :: flxlm_intrp   ! interpolated flux limiter
REAL(KIND=double)                :: beta          ! u/c
REAL(KIND=double)                :: gr_bend2      ! factor to compute gravitational bending of neutrinos from limb
REAL(KIND=double)                :: mu            ! cosine of angle neutrinos from limb with gravitational bending
REAL(KIND=double)                :: mu0           ! fluid frame mu
REAL(KIND=double)                :: flxlm_geom    ! geometrical flux limiter
REAL(KIND=double)                :: flxlm         ! net flux limiter

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: r2     ! r^{2}
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rjmh2  ! rjmh^{2}

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in diffc')
 2001 FORMAT (' Deallocation problem for array ',a10,' in diffc')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (r2(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r2        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rjmh2(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rjmh2     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize r2 and rjmh2
!-----------------------------------------------------------------------

r2                   = zero
rjmh2                = zero

!-----------------------------------------------------------------------
!  Compute r2 and rjmh2
!-----------------------------------------------------------------------

DO j = jr_min,jr_max+1
  r2(j)              = r(j) * r(j)
  rjmh2(j)           = rjmh(j) * rjmh(j)
END DO

!-----------------------------------------------------------------------
!
!         \\\\\ COMPUTE DC WHEN PSI0(J+1,K,N) IS PERTURBED /////
!
!           dc_jmh
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Begin loops
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Interpolated flux limiter
!-----------------------------------------------------------------------

      psi0j          = ( ( rjmh2(j+1) - r2(j) ) * psi0(j  ,k,n) + ( r2(j) - rjmh2(j  ) ) * p_dpsi0(j+1,k,n) ) &
&                    / ( rjmh2(j+1) - rjmh2(j) )
      dpsidr         = ( p_dpsi0(j+1,k,n) - psi0(j,k,n) ) * drjmh_inv(j)
      guiuv          = dpsidr/( psi0j + epsilon )
      tmfp_jj        = DSQRT( DABS( tmfp(j,k,n) * tmfp(j+1,k,n) ) + epsilon )
      flxlm_intrp    = one/( one + third * tmfp_jj * DABS( guiuv ) )

!-----------------------------------------------------------------------
!  Geometrical flux limiter
!-----------------------------------------------------------------------

      flxlm_geom     = one

      IF ( r(j) > r_sphere(k,n,ij_ray,ik_ray)      .and.               &
&          r_sphere(k,n,ij_ray,ik_ray) /= zero            ) THEN
        beta         = u(j)/cvel
        IF ( ireltrns == 0 ) THEN
          gr_bend2   = one
        ELSE
          gr_bend2   = ( one - 2.d0 * g * grvmss(j)/( r(j) * cvel**2 ) ) &
&                    /  ( one - 2.d0 * g * m_sphere(k,n,ij_ray,ik_ray)/( r_sphere(k,n,ij_ray,ik_ray) * cvel**2 ) )
        END IF
        mu           = DSQRT( DMAX1( one - ( r_sphere(k,n,ij_ray,ik_ray)/r(j) )**2 * gr_bend2, epsilon ) )
        mu0          = ( mu - beta )/( one - mu * beta )
        
        flxlm_geom   = half * ( one + mu0 ) * psi0j/( third * tmfp_jj * DABS(dpsidr) + epsilon )

!-----------------------------------------------------------------------
!  Net flux limiter
!-----------------------------------------------------------------------

        flxlm        = DMIN1( flxlm_intrp, flxlm_geom )

      ELSE

        flxlm        = flxlm_intrp

      END IF ! r(j) > r_sphere(k,n,ij_ray,ik_ray), r_sphere(k,n,ij_ray,ik_ray) /= 0

!-----------------------------------------------------------------------
!  Diffusion coefficient
!-----------------------------------------------------------------------

      dc_jph(j,k,n)  = third * tmfp_jj * flxlm

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  End loops
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!
!         \\\\\ COMPUTE DC WHEN PSI0(J  ,K,N) IS PERTURBED /////
!
!           dc_jmh
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Begin loops
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Interpolated flux limiter
!-----------------------------------------------------------------------

      psi0j          = ( ( rjmh2(j+1) - r2(j) ) * p_dpsi0(j  ,k,n) + ( r2(j) - rjmh2(j  ) ) * psi0(j+1,k,n) ) &
&                    / ( rjmh2(j+1) - rjmh2(j) )
      dpsidr         = ( psi0(j+1,k,n) - p_dpsi0(j,k,n) ) * drjmh_inv(j)
      guiuv          = dpsidr/( psi0j + epsilon )
      tmfp_jj        = DSQRT( DABS( tmfp(j,k,n) * tmfp(j+1,k,n) ) + epsilon )
      flxlm_intrp    = one/( one + third * tmfp_jj * DABS( guiuv ) )

!-----------------------------------------------------------------------
!  Geometrical flux limiter
!-----------------------------------------------------------------------

      flxlm_geom     = one

      IF ( r(j) > r_sphere(k,n,ij_ray,ik_ray)  .and.  r_sphere(k,n,ij_ray,ik_ray) /= zero ) THEN
        beta         = u(j)/cvel
        IF ( ireltrns == 0 ) THEN
          gr_bend2   = one
        ELSE
          gr_bend2   = ( one - 2.d0 * g * grvmss(j)/( r(j) * cvel**2 ) ) &
&                    /  ( one - 2.d0 * g * m_sphere(k,n,ij_ray,ik_ray)/( r_sphere(k,n,ij_ray,ik_ray) * cvel**2 ) )
        END IF
        mu           = DSQRT( DMAX1( one - ( r_sphere(k,n,ij_ray,ik_ray)/r(j) )**2 * gr_bend2, epsilon ) )
        mu0          = ( mu - beta )/( one - mu * beta )
        
        flxlm_geom   = half * ( one + mu0 ) * psi0j/( third * tmfp_jj * DABS(dpsidr) + epsilon )

!-----------------------------------------------------------------------
!  Net flux limiter
!-----------------------------------------------------------------------

        flxlm        = DMIN1( flxlm_intrp, flxlm_geom )

      ELSE

        flxlm        = flxlm_intrp

      END IF ! r(j) > r_sphere(k,n,ij_ray, ik_ray), r_sphere(k,n,ij_ray, ik_ray) /= 0

!-----------------------------------------------------------------------
!  Diffusion coefficient
!-----------------------------------------------------------------------

      dc_jmh(j,k,n)  = third * tmfp_jj * flxlm

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  End loops
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    END DO
  END DO
END DO

!........Deallocate arrays..............................................

DEALLOCATE (r2, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r2        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rjmh2, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rjmh2     '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE diffc_d

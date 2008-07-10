SUBROUTINE diffc( jr_min, jr_max, ij_ray, ik_ray, r, u, nx )
!-----------------------------------------------------------------------
!
!    File:         diffc
!    Module:       diffc
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
!  jr_min     : inner radial zone for computing the diffusion coefficient.
!  jr_max     : outer radial zone for computing the diffusion coefficient.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!  r          : radius (cm)
!  u          : x-velocity (cm s^{-1})
!  nx         : x-array extent
!
!    Output arguments:
!  dc(j,k,n)  : diffusion coefficient for neutrinos at j, energy zone k, type n
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, mdl_cnfg_module, nu_dist_module,
!  nu_energy_grid_module prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nnu, ij_ray_dim
USE numerical_module, ONLY : zero, half, one, third, epsilon
USE physcnst_module, ONLY : cvel, g

USE cycle_module
USE edit_module, ONLY : nlog
USE mdl_cnfg_module, ONLY : grvmss
USE nu_dist_module, ONLY : psi0, tmfp, r_sphere, m_sphere, rjmh, drjmh_inv, &
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

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: r      ! radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: u      ! x-velocity (cm s^{-1})

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name
CHARACTER (len=10)               :: outfile       ! character string containing name of model files

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status
INTEGER, PARAMETER               :: iunit=42
INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

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

i_ray                  = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

WRITE (outfile,'(a5,i3.3,a2)') 'diffc', &
&  i_ray,'.d'
OPEN (UNIT=iunit,FILE=outfile,STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=iunit,FILE=outfile,STATUS='old',POSITION='append')

r2                     = zero
rjmh2                  = zero

!-----------------------------------------------------------------------
!  Compute r2 and rjmh2
!-----------------------------------------------------------------------

r2(jr_min:jr_max+1)    = r(jr_min:jr_max+1)    * r(jr_min:jr_max+1)
rjmh2(jr_min:jr_max+1) = rjmh(jr_min:jr_max+1) * rjmh(jr_min:jr_max+1)

DO j = jr_min,jr_max+1
  WRITE (iunit,3001) j, r2(j), rjmh2(j)
 3001 FORMAT (' j=',i3,' r2(j)=',es22.14,' rjmh2(j)=',es22.14)
END DO

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Begin loops
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO j = jr_min,jr_max-1

!-----------------------------------------------------------------------
!  Interpolated flux limiter
!-----------------------------------------------------------------------

!      psi0j            = ( ( rjmh2(j+1) - r2(j) ) * psi0(j  ,k,n) + ( r2(j) - rjmh2(j  ) ) * psi0(j+1,k,n) ) &
!&                      / ( rjmh2(j+1) - rjmh2(j) )
      psi0j            = psi0(j,k,n) * rjmh2(j)/r2(j)
      dpsidr           = ( psi0(j+1,k,n) - psi0(j,k,n) ) * drjmh_inv(j)
      guiuv            = dpsidr/( psi0j + epsilon )
      tmfp_jj          = DSQRT( DABS( tmfp(j,k,n) * tmfp(j+1,k,n) ) + epsilon )
      flxlm_intrp      = one/( one + third * tmfp_jj * DABS( guiuv ) )
  WRITE (iunit,3002) ncycle,j,k,n,psi0j
 3002 FORMAT (' ncycle=',i8,' j=',i3,' k=',i3,' n=',i3,' psi0j=',es22.14)
  WRITE (iunit,3003) j,k,n,dpsidr,psi0(j+1,k,n),psi0(j,k,n),drjmh_inv(j)
 3003 FORMAT (' j=',i3,' k=',i3,' n=',i3,' dpsidr=',es22.14,' psi0(j+1,k,n)=',es22.14, &
 & ' psi0(j,k,n)=',es22.14,' drjmh_inv(j)=',es22.14)
  WRITE (iunit,3004) j,k,n,guiuv
 3004 FORMAT (' j=',i3,' k=',i3,' n=',i3,' guiuv=',es22.14)
  WRITE (iunit,3005) j,k,n,tmfp_jj
 3005 FORMAT (' j=',i3,' k=',i3,' n=',i3,' tmfp_jj=',es22.14)
  WRITE (iunit,3006) j,k,n,flxlm_intrp
 3006 FORMAT (' j=',i3,' k=',i3,' n=',i3,' flxlm_intrp=',es22.14)

!-----------------------------------------------------------------------
!  Geometrical flux limiter
!-----------------------------------------------------------------------

      flxlm_geom       = one

      IF ( r(j) > r_sphere(k,n,ij_ray,ik_ray)  .and.  r_sphere(k,n,ij_ray,ik_ray) /= zero ) THEN
        beta           = u(j)/cvel
        IF ( ireltrns == 0 ) THEN
          gr_bend2     = one
        ELSE
          gr_bend2     = ( one - 2.d0 * g * grvmss(j)/( r(j) * cvel**2 ) ) &
&                      /  ( one - 2.d0 * g * m_sphere(k,n,ij_ray,ik_ray)/( r_sphere(k,n,ij_ray,ik_ray) * cvel**2 ) )
        END IF
        mu             = DSQRT( DMAX1( one - ( r_sphere(k,n,ij_ray,ik_ray)/r(j) )**2 * gr_bend2, epsilon ) )
        mu0            = ( mu - beta )/( one - mu * beta )
  WRITE (iunit,3007) j,k,n,mu,u(j)
 3007 FORMAT (' j=',i3,' k=',i3,' n=',i3,' mu=',es22.14,' u(j)=',es22.14)
  WRITE (iunit,3008) j,k,n,mu0
 3008 FORMAT (' j=',i3,' k=',i3,' n=',i3,' mu0=',es22.14)
        
        flxlm_geom     = half * ( one + mu0 ) * psi0j/( third * tmfp_jj * DABS(dpsidr) + epsilon )
  WRITE (iunit,3009) j,k,n,flxlm_geom, half * ( one + mu0 )/( third * tmfp_jj + epsilon )
 3009 FORMAT (' j=',i3,' k=',i3,' n=',i3,' flxlm_geom=',es22.14,' half * ( one + mu0 )/( third * tmfp_jj + epsilon )=',es22.14)

!-----------------------------------------------------------------------
!  Net flux limiter
!-----------------------------------------------------------------------

        flxlm          = DMIN1( flxlm_intrp, flxlm_geom )

      ELSE

        flxlm          = flxlm_intrp

      END IF ! r(j) > r_sphere(k,n,ij_ray,ik_ray), r_sphere(k,n,ij_ray,ik_ray) /= 0

!-----------------------------------------------------------------------
!  Diffusion coefficient
!-----------------------------------------------------------------------

      dc(j,k,n)        = third * tmfp_jj * flxlm
  WRITE (iunit,3010) j,k,n,flxlm_intrp,flxlm_geom,flxlm,tmfp_jj,dc(j,k,n)
 3010 FORMAT (' j=',i3,' k=',i3,' n=',i3,' flxlm_intrp=',es22.14,' flxlm_geom=',es22.14, &
 & ' flxlm=',es22.14,' tmfp_jj=',es22.14,' dc(j,k,n)=',es22.14)

!-----------------------------------------------------------------------
!  Store flux limiter
!-----------------------------------------------------------------------

      fluxlm(j,k,n)    = flxlm
 
!-----------------------------------------------------------------------
!  Store zone-edged mean free path
!-----------------------------------------------------------------------

      tmfp_j(j,k,n)    = tmfp_jj

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  End loops
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    END DO
  END DO
END DO

dc    (jr_max,:,:)     = dc    (jr_max-1,:,:)
fluxlm(jr_max,:,:)     = fluxlm(jr_max-1,:,:)

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (r2, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r2        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rjmh2, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rjmh2     '; WRITE (nlog,2001) var_name; END IF

CLOSE (UNIT=iunit,status='keep')

RETURN
END SUBROUTINE diffc

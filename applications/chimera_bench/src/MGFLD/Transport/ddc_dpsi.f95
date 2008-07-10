SUBROUTINE ddc_dpsi( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& r, ddcpsjph, ddcpsjmh, nx, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         ddc_dpsi
!    Module:       ddc_dpsi
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/12/00
!
!    Purpose:
!      To compute derivatives of the neutrino diffusion coefficients with
!       respect to the neutrino occupation functions.
!
!    Subprograms called:
!  diff         : computes diffusion coefficients for the current time step
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
!  r            : radius (cm)
!  nx           : x-array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!
!    Output arguments:
!  ddcpsjph     : d(dc)/d(psi(j+1,k,n)
!  ddcpsjmh     : d(dc)/d(psi(j  ,k,n)
!
!    Input arguments (common):
!      none
!
!    Output arguments (common):
!  psi1(j,k,n)  : first angular moment of theneutrino distribution function
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, nu_dist_module, nu_energy_grid_module, prb_cntl_module,
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, epsilon

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : psi0
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

REAL(KIND=double), INTENT(in), DIMENSION(nx)          :: r         ! radius (cm)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu) :: ddcpsjph  ! derivative of dc wrt psi0(j+1,k,n)
REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu) :: ddcpsjmh  ! derivative of dc wrt psi0(j  ,k,n)

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: p_pdpsi0     !
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: p_mdpsi0     !
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dc_pjph      !
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dc_pjmh      !
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dc_mjph      !
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dc_mjmh      !

 1001 FORMAT (' Allocation problem for array ',a10,' in ddc_dpsi')
 2001 FORMAT (' Deallocation problem for array ',a10,' in ddc_dpsi')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( idiff /= 0 ) RETURN

!........Allocate arrays................................................

ALLOCATE (p_pdpsi0(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_pdpsi0  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (p_mdpsi0(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_mdpsi0  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dc_pjph(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_pjph   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dc_pjmh(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_pjmh   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dc_mjph(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_mjph   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dc_mjmh(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_mjmh   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

p_pdpsi0                = zero
p_mdpsi0                = zero
dc_pjph                 = zero
dc_pjmh                 = zero
dc_mjph                 = zero
dc_mjmh                 = zero


!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE DC WHEN PSI0 IS INCREMENTED /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Increment the neutrino occupation function
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      DO j = jr_min,jr_max
        p_pdpsi0(j,k,n) = psi0(j,k,n) + DMAX1( 0.05d0 * psi0(j,k,n), 1.d-06 )
      END DO ! j
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!  Compute dc_pjph and dc_pjmh
!-----------------------------------------------------------------------

CALL diffc_d( jr_min, jr_max, ij_ray, ik_ray, r, psi0, p_pdpsi0, dc_pjph, &
& dc_pjmh, nx, nez, nnu )

!-----------------------------------------------------------------------
!
!          \\\\\ COMPUTE DC WHEN PSI0 IS DECREMENTED /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Decrement the neutrino occupation function
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      DO j = jr_min,jr_max
        p_mdpsi0(j,k,n) = psi0(j,k,n) - DMAX1( 0.05d0 * psi0(j,k,n), 1.d-06 )
      END DO ! j
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n

!-----------------------------------------------------------------------
!  Compute dc_mjph and dc_mjmh
!-----------------------------------------------------------------------

CALL diffc_d( jr_min, jr_max, ij_ray, ik_ray, r, psi0, p_mdpsi0, dc_mjph, &
& dc_mjmh, nx, nez, nnu )

!-----------------------------------------------------------------------
!
!             \\\\\ COMPUTE DERIVATIVES OF THE DC /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute dc_pjph and ddcpsjmh
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) /= 0 ) THEN
    DO k = 1,nnugp(n)
      DO j = jr_min,jr_max
        ddcpsjph(j,k,n) = + 0.d0 * ( dc_pjph(j,k,n) - dc_mjph(j,k,n) ) &
&                       / ( p_pdpsi0(j+1,k,n) - p_mdpsi0(j+1,k,n) + epsilon )
        ddcpsjmh(j,k,n) = + 0.d0 * ( dc_pjmh(j,k,n) - dc_mjmh(j,k,n) ) &
&                       / ( p_pdpsi0(j  ,k,n) - p_mdpsi0(j  ,k,n) + epsilon )
      END DO ! j
    END DO ! k
  END IF ! nnugp(n) /= 0
END DO ! n
  WRITE (nlog,3001) (dc_pjph(jr_min,k,1),k=1,20)
 3001 FORMAT (' dc_pjph(jr_min,k,1)=',20es11.3)
  WRITE (nlog,3002) (dc_mjph(jr_min,k,1),k=1,20)
 3002 FORMAT (' dc_mjph(jr_min,k,1)=',20es11.3)
  WRITE (nlog,3003) (p_pdpsi0(jr_min,k,1),k=1,20)
 3003 FORMAT (' p_pdpsi0(jr_min,k,1)=',20es11.3)
  WRITE (nlog,3004) (p_mdpsi0(jr_min,k,1),k=1,20)
 3004 FORMAT (' p_mdpsi0(jr_min,k,1)=',20es11.3)

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (p_pdpsi0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_pdpsi0  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (p_mdpsi0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'p_mdpsi0  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dc_pjph, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_pjph   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dc_pjmh, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_pjmh   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dc_mjph, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_mjph   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dc_mjmh, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dc_mjmh   '; WRITE (nlog,2001) var_name; END IF


RETURN
END SUBROUTINE ddc_dpsi

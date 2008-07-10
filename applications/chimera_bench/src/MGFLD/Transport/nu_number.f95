SUBROUTINE nu_number( jr_min, jr_max, n, ij_ray, ik_ray, nx, nez, nnu, r,&
& u, psi0, psi1 )
!-----------------------------------------------------------------------
!
!    File:         nu_number
!    Module:       nu_number
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/25/01
!
!    Purpose:
!      To compute nnujcr(j,n), unujcr(j,n), nnucr(n,ij_ray,ik_ray), and unucr(n).
!       The neutrino number and energies are computed assuming a spherically
!       symmetric shell. They should be multiplied by the solid angle
!       subtended by a radial ray to obtain the actual number and energy.
!
!    Variables that must be passed through common:
!  ncoefa(j,k)    : 4.*pi/((h*c)**3)*w**2*dw
!  ecoefa(j,k)    : 4.*pi/((h*c)**3)*w**3*dw
!
!    Subprograms called:
!        pre_trans
!
!    Input arguments:
!  jr_min         : inner radial zone of region for which configuration edit is to be made.
!  jr_max         : outer radial zone of region for which configuration edit is to be made.
!  n              : neutrino flavor index
!  ij_ray         : j-index of a radial ray
!  ik_ray         : k-index of a radial ray
!  nx             : x-array extent
!  nez            : neutrino energy array extent
!  nnu            : neutrino flavor array extent
!  r              : radius (cm)
!  u              : velocity x-components (cm s^{-1})
!  psi0           : zero angular moments of the neutrino occupation number
!  psi1           : first angular moments of the neutrino occupation number
!
!    Output arguments:
!        none
!
!    Output arguments (common):
!  nnujcr(j,n)    : net number of n-type neutrinos currently residing in radial zone j
!  unujcr(j,n)    : total energy of n-type neutrinos currently residing in radial zone j (ergs)
!  unujinfty(j,n) : total energy of n-type neutrinos currently residing in radial zone j as seen from infity(ergs)
!  nnucr(n,ij_ray,ik_ray)
!                 : net number of n-type neutrinos currently residing in the core
!  unucr(n)       : total energy of n-type neutrinos currently residing in the core (ergs)
!  unuinfty(n)    : total energy of n-type neutrinos currently residing in the core as seen from infinity (ergs)
!
!    Modules used:
!  kind_module, numerical_module, physcnst_module
!  nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, half, frpith
USE physcnst_module, ONLY: cvel

USE nu_dist_module, ONLY: ncoefa, ecoefa, stwt, unujcr, unucr, unujinfty, unuinfty, &
&   nnujcr, nnucr, gamgr_nu
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: n             ! neutrino type index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: r    ! radius (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: u    ! zone-centered x-velocity velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu) :: psi0 ! zero moment of the neutrino occupation probability
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu) :: psi1 ! first moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! energy zone index

REAL(KIND=double)                :: rnnutt        ! temporary variable for computing neutrino number
REAL(KIND=double)                :: unucrtt       ! temporary variable for computing neutrino energy
REAL(KIND=double)                :: psi0_infty    ! psi0 as observed from infinity

REAL(KIND=double), DIMENSION(nx) :: vol           ! radius (cm)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                     \\\\\ INITIALIZE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute volumes
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  vol(j)                   = frpith * ( r(j) - r(j-1) )                      &
&                          * ( r(j) * ( r(j) - r(j-1) ) + r(j-1) * r(j-1) )
END DO

!-----------------------------------------------------------------------
!  Initialize arrays
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  unujcr(j,n)              = zero
  unujinfty(j,n)           = zero
  nnujcr(j,n)              = zero
END DO

unucr(n)                   = zero
unuinfty(n)                = zero
nnucr(n,ij_ray,ik_ray)     = zero

!-----------------------------------------------------------------------
!
!                   \\\\\ UNUJCR AND UNUCR /////
!
!-----------------------------------------------------------------------

DO k = 1,nnugp(n)
  DO j = jr_min,jr_max
    unucrtt                = ecoefa(j,k) * psi0(j,k,n) * stwt(n) * vol(j)
    unujcr(j,n)            = unujcr(j,n) + unucrtt
    unucr(n)               = unucr(n) + unucrtt
  END DO
END DO

!-----------------------------------------------------------------------
!
!                \\\\\ UNUJINFTY AND UNUINFTY /////
!
!-----------------------------------------------------------------------

DO k = 1,nnugp(n)
  DO j = jr_min,jr_max
    psi0_infty             = half * ( gamgr_nu(j) + gamgr_nu(j-1) ) * psi0(j,k,n) &
&                          + ( half * ( u(j) + u(j-1) )/cvel ) * psi1(j,k,n)
    unucrtt                = ecoefa(j,k) * stwt(n) * vol(j) * psi0_infty
    unujinfty(j,n)         = unujinfty(j,n) + unucrtt
    unuinfty(n)            = unuinfty(n) + unucrtt
  END DO
END DO

!-----------------------------------------------------------------------
!
!                  \\\\\ NNUJCR AND NNUCR /////
!
!-----------------------------------------------------------------------

DO k = 1,nnugp(n)
  DO j = jr_min,jr_max
    rnnutt                 = ncoefa(j,k) * psi0(j,k,n) * stwt(n) * vol(j)
    nnujcr(j,n)            = nnujcr(j,n) + rnnutt
    nnucr(n,ij_ray,ik_ray) = nnucr(n,ij_ray,ik_ray) + rnnutt
  END DO
END DO

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE nu_number

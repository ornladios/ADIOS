SUBROUTINE pblmst1( ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         pblmst1
!    Module:       pblmst1
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/29/05
!
!    Purpose:
!        To perturb the density near the polar axis to cause an asymmetric
!         collapse.
!
!    Input arguments:
!  ij_ray : j-index of a radial ray
!  ik_ray : k-index of a radial ray
!
!    Include files:
!  cycle_module, mdl_cnfg_module, parallel_module, radial_ray_module
!     
!-----------------------------------------------------------------------

USE cycle_module, ONLY : ncycle
USE mdl_cnfg_module, ONLY : rhor, rho
USE parallel_module, ONLY : myid
USE radial_ray_module, ONLY : imin, imax, rho_c, rho_ci

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i               ! do index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( ncycle == 0 ) THEN
  IF ( myid <= 4 ) THEN
    DO i = imin,imax
      rho_c (i,ij_ray,ik_ray) = 1.05d0 * rho_c (i,ij_ray,ik_ray)
      rho_ci(i,ij_ray,ik_ray) = 1.05d0 * rho_ci(i,ij_ray,ik_ray)
      rhor  (i+1)             = 1.05d0 * rho_c (i,ij_ray,ik_ray)
      rho   (i+1)             = 1.05d0 * rho_c (i,ij_ray,ik_ray)
    END DO ! i
  END IF ! myid <= 4
END IF ! ncycle == 0

RETURN
END SUBROUTINE pblmst1

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
!        To set the velocities to zero for performing a static test.
!
!    Input arguments:
!  ij_ray : j-index of a radial ray
!  ik_ray : k-index of a radial ray
!
!    Include files:
!  numerical_module
!  mdl_cnfg_module, radial_ray_module
!     
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE mdl_cnfg_module, ONLY : u
USE radial_ray_module, ONLY : imin, imax, u_c, u_ci

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

DO i = imin,imax
  u_c (i,ij_ray,ik_ray) = zero
  u_ci(i,ij_ray,ik_ray) = zero
  u   (i+1)             = zero
END DO ! i

RETURN
END SUBROUTINE pblmst1

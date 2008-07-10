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
!  array_module
!  cycle_module, mdl_cnfg_module, parallel_module, radial_ray_module
!     
!-----------------------------------------------------------------------

USE array_module, ONLY : n_proc_y, ij_ray_dim, ik_ray_dim

USE cycle_module, ONLY : ncycle
USE mdl_cnfg_module, ONLY : rhor, rho, t, tr, ye, yer, u
USE parallel_module, ONLY : myid
USE radial_ray_module, ONLY : imin, imax, rho_c, rho_ci, t_c, t_ci, ye_c, &
& ye_ci, u_c, u_ci

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

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

j_ray              = MOD( myid, n_proc_y ) * ij_ray_dim + ij_ray
k_ray              = ( myid/n_proc_y ) * ik_ray_dim + ik_ray

IF ( ncycle == 0 ) THEN
  IF ( k_ray == 2 ) THEN
    DO i = imin,imax
      rho_c (i,ij_ray,ik_ray) = 1.05d0 * rho_c (i,ij_ray,ik_ray)
      rho_ci(i,ij_ray,ik_ray) = 1.05d0 * rho_c (i,ij_ray,ik_ray)
      rhor  (i+1)             = 1.05d0 * rho_c (i,ij_ray,ik_ray)
      rho   (i+1)             = 1.05d0 * rho_c (i,ij_ray,ik_ray)
      t_c   (i,ij_ray,ik_ray) = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      t_ci  (i,ij_ray,ik_ray) = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      tr    (i+1)             = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      t     (i+1)             = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      t_c   (i,ij_ray,ik_ray) = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      t_ci  (i,ij_ray,ik_ray) = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      tr    (i+1)             = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      t     (i+1)             = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      t     (i+1)             = 1.05d0 * t_c   (i,ij_ray,ik_ray)
      ye_c  (i,ij_ray,ik_ray) = 1.05d0 * ye_c  (i,ij_ray,ik_ray)
      ye_ci (i,ij_ray,ik_ray) = 1.05d0 * ye_c  (i,ij_ray,ik_ray)
      yer   (i+1)             = 1.05d0 * ye_c  (i,ij_ray,ik_ray)
      ye    (i+1)             = 1.05d0 * ye_c  (i,ij_ray,ik_ray)
      u_c   (i,ij_ray,ik_ray) = 1.05d0 * u_c   (i,ij_ray,ik_ray)
      u_ci  (i,ij_ray,ik_ray) = 1.05d0 * u_c   (i,ij_ray,ik_ray)
      u     (i+1)             = 1.05d0 * u_c   (i,ij_ray,ik_ray)
      v_c   (i,ij_ray,ik_ray) = 1.d+03
      v_ci  (i,ij_ray,ik_ray) = 1.d+03
      v     (i+1)             = 1.d+03
      w_c   (i,ij_ray,ik_ray) = 1.d+03
      w_ci  (i,ij_ray,ik_ray) = 1.d+03
      w     (i+1)             = 1.d+03
    END DO ! i
  END IF ! myid <= 4
END IF ! ncycle == 0

RETURN
END SUBROUTINE pblmst1

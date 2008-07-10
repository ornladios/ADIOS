SUBROUTINE ye_fix( nx, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         ye_fix
!    Module:       ye_fix
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/31/05
!
!    Purpose:
!      To set the value of ye for nonNSE material to values consistent
!       the chemical composition.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x-array array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  radial_ray_module, 
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : epsilon

USE radial_ray_module, ONLY : imin, imax, nuc_number, nse_c, ye_c, ye_ci, &
& xn_c, a_nuc_c, z_nuc_c, z_nuc_rep_c, a_nuc_rep_c

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx           ! energy array dimension
INTEGER, INTENT(in)              :: ij_ray_dim   ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim   ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ij_ray       ! j-index of a radial ray
INTEGER                          :: ik_ray       ! k-index of a radial ray
INTEGER                          :: i            ! radial zone index
INTEGER                          :: n_nucp1      ! radial zone index

REAL(KIND=double)                :: z_tot        ! total proton number (used for computing ye)
REAL(KIND=double)                :: a_tot        ! total mass number (used for computing ye)
REAL(KIND=double)                :: ye_tot       ! electron fraction

!-----------------------------------------------------------------------
!
!              \\\\\ SET NES FALG FOR OUTER GHOSTS /////
!
!-----------------------------------------------------------------------

n_nucp1                        = nuc_number + 1

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    DO i = imin,imax
      IF ( nse_c(i,ij_ray,ik_ray) == 0 ) THEN
        z_tot                  = SUM( xn_c(i,nuc_number,ij_ray,ik_ray) * z_nuc_c(1:nuc_number) )
        a_tot                  = SUM( xn_c(i,nuc_number,ij_ray,ik_ray) * a_nuc_c(1:nuc_number) )
        z_tot                  = z_tot + xn_c(i,n_nucp1,ij_ray,ik_ray) * z_nuc_rep_c(i,ij_ray,ik_ray)
        a_tot                  = a_tot + xn_c(i,n_nucp1,ij_ray,ik_ray) * a_nuc_rep_c(i,ij_ray,ik_ray)
        ye_tot                 = z_tot/( a_tot + epsilon )
        ye_c(i,ij_ray,ik_ray)  = ye_tot
        ye_ci(i,ij_ray,ik_ray) = ye_tot
      END IF ! nse_c(imax,ij_ray,ik_ray) == 0
    END DO ! i = imin,imax
  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

RETURN
END SUBROUTINE ye_fix

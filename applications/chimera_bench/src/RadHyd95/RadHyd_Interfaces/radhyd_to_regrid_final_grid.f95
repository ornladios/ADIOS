SUBROUTINE radhyd_to_regrid_final_grid( nx, ij_ray_dim, ny, ik_ray_dim, &
& nz )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_regrid_final_grid
!    Module:       radhyd_to_regrid_final_grid
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/12/07
!
!    Purpose:
!      To port radial_ray_module variables into the final grid regridder
!       subroutine, and regridded variabbles back to the radial_ray_module.
!
!    Input arguments:
!  nx                : x-array extent
!  ij_ray_dim        : the number of y-zones on a processor before swapping
!                       with y
!  ny                : y-array extent
!  ik_ray_dim        : the number of z-zones on a processor before swapping
!                       with z
!  nz                : z-array extent
!
!    Subprograms called:
!  regrid_final_grid : executes the remap along the x-direction
!
!    Include files:
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, nse_c, x_ef, dx_cf, x_cf, rhobar, &
& regrid, rho_regrid, grid_frac, int_pre_b, int_post_b, time, t_bounce

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! the number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: ik_ray_dim    ! the number of z-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: nz            ! z-array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer variables to remap arrays
!-----------------------------------------------------------------------

CALL regrid_final_grid( imin, imax, nx, ij_ray_dim, ny, ik_ray_dim, nz, &
& nse_c, x_ef, dx_cf, x_cf, rhobar, regrid, rho_regrid, grid_frac, &
& int_pre_b, int_post_b, time, t_bounce )

RETURN
END SUBROUTINE radhyd_to_regrid_final_grid

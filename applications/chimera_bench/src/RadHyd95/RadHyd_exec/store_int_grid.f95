 SUBROUTINE store_int_grid
!-----------------------------------------------------------------------
!
!    File:         store_int_grid
!    Module:       store_int_grid
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/05/05
!
!    Purpose:
!      To store initial grid before the Lagrangian hydro step.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx  : x-array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : x_ef, dx_cf, x_cf, y_ef, dy_cf, y_cf, z_ef, &
& dz_cf, z_cf, x_ei, dx_ci, x_ci, y_ei, dy_ci, y_ci, z_ei, dz_ci, z_ci
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Store initial grid
!-----------------------------------------------------------------------

x_ei               = x_ef
y_ei               = y_ef
z_ei               = z_ef

dx_ci              = dx_cf
dy_ci              = dy_cf
dz_ci              = dz_cf

x_ci               = x_cf
y_ci               = y_cf
z_ci               = z_cf

RETURN
END SUBROUTINE store_int_grid

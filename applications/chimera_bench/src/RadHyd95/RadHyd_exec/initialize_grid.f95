SUBROUTINE initialize_grid
!-----------------------------------------------------------------------
!
!    File:         initialize_grid
!    Module:       initialize_grid
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/05/05
!
!    Purpose:
!      To initial grid before the probklem execution.
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
!  Initialize final grid to initial grid
!-----------------------------------------------------------------------

x_ef               = x_ei
y_ef               = y_ei
z_ef               = z_ei

dx_cf              = dx_ci
dy_cf              = dy_ci
dz_cf              = dz_ci

x_cf               = x_ci
y_cf               = y_ci
z_cf               = z_ci

RETURN
END SUBROUTINE initialize_grid

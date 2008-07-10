 SUBROUTINE set_final_azimuthal_grid
!-----------------------------------------------------------------------
!
!    File:         set_final_azimuthal_grid
!    Module:       set_final_azimuthal_grid
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/05/05
!
!    Purpose:
!      To Set the angular grid at the end of the hydro step.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  i_ray_dim : number of radial rays on a processor
!  nx        : x-array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : z_ei, dz_ci, z_ci, z_ef, dz_cf, z_cf

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                         \\\\\ EULERIAN /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set final grid to the initial grid
!-----------------------------------------------------------------------

z_ef                      = z_ei
dz_cf                     = dz_ci
z_cf                      = z_ci

RETURN
END SUBROUTINE set_final_azimuthal_grid

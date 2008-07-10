SUBROUTINE radhyd_to_poisson( mode, nx, ij_ray_dim, ik_ray_dim, ny, nz )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_poisson
!    Module:       radhyd_to_poisson
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the x-array Lagrangian hydro modules via subroutine
!       evh1_x_lagr_inout.
!
!    Input arguments:
!  mode       : Legendre polynomial evaluation flag
!  nx         : x-array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  ny         : y-array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  poisson   : computes the gravitational potential and accelerations
!
!    Include files:
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, x_ei, x_ci, dx_ci, y_ei, y_ci, &
& dy_ci, z_ei, dz_ci, rho_c, grav_x_c, grav_y_c, grav_z_c, grav_pot_c,   &
& grav_x_e, grav_y_e, grav_z_e, grav_pot_e, i_grav

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: mode          ! Legendre polynomial evaluation flag
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if i_grav /= 2
!-----------------------------------------------------------------------

IF ( i_grav /= 2 ) RETURN

!-----------------------------------------------------------------------
!  Call Possion solver
!-----------------------------------------------------------------------


CALL poisson( mode, imin, imax, nx, ij_ray_dim, ik_ray_dim, ny, nz,   & 
& x_ei, x_ci, dx_ci, y_ei, y_ci, dy_ci, z_ei, dz_ci, rho_c, grav_x_c, &
& grav_y_c, grav_z_c, grav_pot_c, grav_x_e, grav_y_e, grav_z_e, grav_pot_e )

RETURN
END SUBROUTINE radhyd_to_poisson

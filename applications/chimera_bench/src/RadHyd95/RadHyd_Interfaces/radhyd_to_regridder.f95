SUBROUTINE radhyd_to_regridder( nx, ij_ray_dim, ny, ik_ray_dim, nz, nez, &
& nnu, nnc, l_regrid )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_regridder
!    Module:       radhyd_to_regridder
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/04/07
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the x-array remap modules via subroutine remap_x_inout.
!
!    Input arguments:
!  nx            : x-array extent
!  ij_ray_dim    : the number of y-zones on a processor before swapping with y
!  ny            : y-array extent
!  ik_ray_dim    : the number of z-zones on a processor before swapping with z
!  nz            : z-array extent
!  nez           : neutrino energy array extent
!  nnu           : neutrino flavor array extent
!  nnc           : composition array extent
!
!    Output arguments:
!  l_regrid      : regrid flag
!
!    Subprograms called:
!  remap_x_inout : executes the remap along the x-direction
!
!    Include files:
!  radial_ray_module, evh1_global
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, nse_c, x_ef, dx_cf, x_cf, rhobar, &
& rho_c, t_c, ye_c, u_c, u_e, v_c, w_c, psi0_c, regrid, &
& rho_regrid, time, t_bounce

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
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER                          :: nnc           ! composition array extent

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

LOGICAL, INTENT(out)             :: l_regrid      ! regrid flag

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer variables to remap arrays
!-----------------------------------------------------------------------

CALL regridder( imin, imax, nx, ij_ray_dim, ny, ik_ray_dim, nz, nez, &
& nnu, nnc, nse_c, x_ef, dx_cf, x_cf, rhobar, rho_c, t_c, ye_c, u_c, u_e, &
& v_c, w_c, psi0_c, regrid, rho_regrid, time, t_bounce, l_regrid )

RETURN
END SUBROUTINE radhyd_to_regridder

SUBROUTINE radhyd_to_nu_e_advct_y( nx, ny, ji_ray, jk_ray, j_ray_dim, &
& ik_ray_dim, i_radial, j_radial, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nu_e_advct_y
!    Module:       radhyd_to_nu_e_advct_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the neutrino energy advection modules via subroutine
!       nu_energy_advct_inout_y.
!
!    Subprograms called:
!  nu_energy_advct_inout_y
!
!    Input arguments:
!  nx         : x-array extent
!  ny         : y-array extent
!  ji_ray     : x (radial) index of a specific y (angular) ray
!  jk_ray     : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim  : the number of radial zones on a processor after swapping with y
!  ik_ray_dim : the number of z-zones on a processor before swapping with z
!  i_radial   : the unshifted radial zone corresponding to ji_ray, jk_ray
!  j_radial   : the shifted radial zone corresponding to ji_ray, jk_ray
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  radial_ray_module, angular_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, nprint, ncycle, &
& dtime=>dtnph, time, x_ef, x_cf, y_ei, y_el, rhobar
USE angular_ray_module, ONLY : rho_yi, rho_y, t_y, ye_y, v_y, psi0_y, &
& psi1_y

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ji_ray, jk_ray
INTEGER, INTENT(in)              :: j_radial      ! the shifted radial zone corresponding to ji_ray, jk_ray
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!  Transfer variables to nu_energy_advct_y modules
!-----------------------------------------------------------------------

CALL nu_energy_advct_inout_y( imin, imax, jmin, jmax, nx, ny, ji_ray, jk_ray, &
& j_ray_dim, ik_ray_dim, i_radial, j_radial, nez, nnu, nprint, rho_yi, rho_y, &
& x_ef, x_cf, y_ei, y_el, t_y, ye_y, v_y, psi0_y, psi1_y, dtime, rhobar )

RETURN
END SUBROUTINE radhyd_to_nu_e_advct_y

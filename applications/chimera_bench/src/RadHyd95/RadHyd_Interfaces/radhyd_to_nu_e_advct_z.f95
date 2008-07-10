SUBROUTINE radhyd_to_nu_e_advct_z( nx, nz, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, i_radial, j_radial, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nu_e_advct_z
!    Module:       radhyd_to_nu_e_advct_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/07
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the neutrino energy advection modules via subroutine
!       nu_energy_advct_inout_z.
!
!    Subprograms called:
!  nu_energy_advct_inout_z
!
!    Input arguments:
!  nx         : x-array extent
!  nz         : z-array extent
!  ki_ray     : x (radial) index of a specific z (azimuthal) ray
!  kj_ray     : y (angular) index of a specific z (azimuthal) ray
!  j_ray_dim  : the number of radial zones on a processor after swapping with y
!  ik_ray_dim : the number of z-zones on a processor before swapping with z
!  i_radial   : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  j_radial   : the shifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  radial_ray_module, azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, nprint, &
& ncycle, dtime=>dtnph, time, x_ef, x_cf, z_ei, z_el, rhobar
USE azimuthal_ray_module, ONLY : rho_zi, rho_z, t_z, ye_z, v_z, psi0_z, &
& psi1_z

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: k_ray_dim     ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! unshifted radial zone index for a given ki_ray, kj_ray
INTEGER, INTENT(in)              :: j_radial      ! shifted radial zone index for a given ki_ray, kj_ray
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!  Transfer variables to nu_energy_advct_y modules
!-----------------------------------------------------------------------

CALL nu_energy_advct_inout_z( imin, imax, kmin, kmax, nx, nz, ki_ray, kj_ray, &
& ij_ray_dim, k_ray_dim, i_radial, j_radial, nez, nnu, nprint, rho_zi, rho_z, &
& x_ef, x_cf, z_ei, z_el, t_z, ye_z, v_z, psi0_z, psi1_z, dtime, rhobar )

RETURN
END SUBROUTINE radhyd_to_nu_e_advct_z

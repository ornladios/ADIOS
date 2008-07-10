 SUBROUTINE store_int_azimuthal_var( ki_ray, kj_ray, nz )
!-----------------------------------------------------------------------
!
!    File:         store_int_azimuthal_var
!    Module:       store_int_azimuthal_var
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/23/07
!
!    Purpose:
!      To store coordinate and state variables before the azimuthal
!       Lagrangian hydro step.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ki_ray : x (radial) index of a specific z (azimuthal) ray
!  kj_ray : y (angular) index of a specific z (azimuthal) ray
!  nz     : z-array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE azimuthal_ray_module, ONLY : rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z, &
& rho_zi, t_zi, ye_zi, ei_zi, u_zi, v_zi, w_zi

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: nz              ! z-array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!         \\\\\ STORE INITIAL VARIABLES FOR THE Y-SWEEP /////
!
!-----------------------------------------------------------------------

rho_zi(1:nz,kj_ray,ki_ray) = rho_z(1:nz,kj_ray,ki_ray)
t_zi  (1:nz,kj_ray,ki_ray) = t_z  (1:nz,kj_ray,ki_ray)
ye_zi (1:nz,kj_ray,ki_ray) = ye_z (1:nz,kj_ray,ki_ray)
ei_zi (1:nz,kj_ray,ki_ray) = ei_z (1:nz,kj_ray,ki_ray)
u_zi  (1:nz,kj_ray,ki_ray) = u_z  (1:nz,kj_ray,ki_ray)
v_zi  (1:nz,kj_ray,ki_ray) = v_z  (1:nz,kj_ray,ki_ray)
w_zi  (1:nz,kj_ray,ki_ray) = w_z  (1:nz,kj_ray,ki_ray)

RETURN
END SUBROUTINE store_int_azimuthal_var

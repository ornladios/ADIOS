SUBROUTINE radhyd_to_hydro_z_time_step( kmin, kmax, ki_ray, kj_ray, &
& ij_ray_dim, k_ray_dim, dtime_z_hydro )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_hydro_z_time_step
!    Module:       radhyd_to_hydro_z_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/23/07
!
!    Purpose:
!      To determine the minimum time step for the y-hydro sweep for a
!       given j_ray.
!
!    Input arguments:
!  kmin          : inner z-array index
!  kmax          : outer z-array index
!  ki_ray        : x (radial) index of a specific z (azimuthal) ray
!  kj_ray        : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim    : number of y (angular) zones on a processor before swapping with y
!  k_ray_dim     : number of x (radial) zones on a processor after swapping with z
!
!    Output arguments:
!  dtime_z_hydro : minimum time step by z-hydro criteria
!
!    Subprograms called:
!  courant_z_time_step : selects the minimum time step as given by the
!   Caurant condition applied to the y-hydro sweeo for a given j_ray.
!
!    Include files:
!  kind_module
!  azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE azimuthal_ray_module, ONLY : dt_z, jdt_z

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)            :: kmin          ! inner z-array index
INTEGER, INTENT(in)            :: kmax          ! outer z-array index
INTEGER, INTENT(in)            :: ki_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)            :: kj_ray        ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)            :: ij_ray_dim    ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)            :: k_ray_dim     ! number of radial zones on a processor before swapping with z

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out) :: dtime_z_hydro ! minimum time step by y-hydro criteria

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Deterimine the y-hydro time step for a given j_ray
!-----------------------------------------------------------------------

CALL courant_z_time_step( kmin, kmax, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, dt_z, jdt_z, dtime_z_hydro )

RETURN
END SUBROUTINE radhyd_to_hydro_z_time_step

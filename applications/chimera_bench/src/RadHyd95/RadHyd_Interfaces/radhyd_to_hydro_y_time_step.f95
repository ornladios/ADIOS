SUBROUTINE radhyd_to_hydro_y_time_step( jmin, jmax, ji_ray, jk_ray, &
& j_ray_dim, ik_ray_dim, dtime_y_hydro )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_hydro_y_time_step
!    Module:       radhyd_to_hydro_y_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/05/05
!
!    Purpose:
!      To determine the minimum time step for the y-hydro sweep for a
!       given j_ray.
!
!    Input arguments:
!  jmin          : inner y-array index
!  jmax          : outer y-array index
!  ji_ray        : x (radial) index of a specific y (angular) ray
!  jk_ray        : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim     : the number of radial zones on a processor after swapping with y
!  ik_ray_dim    : the number of z-zones on a processor before swapping with z
!
!    Output arguments:
!  dtime_y_hydro : minimum time step by y-hydro criteria
!
!    Subprograms called:
!  courant_y_time_step : selects the minimum time step as given by the
!   Caurant condition applied to the y-hydro sweeo for a given j_ray.
!
!    Include files:
!  kind_module
!  angular_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE angular_ray_module, ONLY : dt_y, jdt_y

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)            :: jmin          ! inner y-array index
INTEGER, INTENT(in)            :: jmax          ! outer y-array index
INTEGER, INTENT(in)            :: ji_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)            :: jk_ray        ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)            :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)            :: ik_ray_dim    ! number of radial zones on a processor before swapping with z

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out) :: dtime_y_hydro ! minimum time step by y-hydro criteria

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Deterimine the y-hydro time step for a given j_ray
!-----------------------------------------------------------------------

CALL courant_y_time_step( jmin, jmax, ji_ray, jk_ray, j_ray_dim, &
& ik_ray_dim, dt_y, jdt_y, dtime_y_hydro )

RETURN
END SUBROUTINE radhyd_to_hydro_y_time_step

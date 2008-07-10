SUBROUTINE courant_y_time_step( jmin, jmax, ji_ray, jk_ray, j_ray_dim, &
& ik_ray_dim, dt_y, jdt_y, dtime_y_hydro )
!-----------------------------------------------------------------------
!
!    File:         courant_y_time_step
!    Module:       courant_y_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To select the time step restricted by matter processes during
!       the x-sweep.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jmin          : inner y-array index
!  jmax          : outer y-array index
!  ji_ray        : x (radial) index of a specific y (angular) ray
!  jk_ray        : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim     : the number of radial zones on a processor after swapping with y
!  ik_ray_dim    : the number of z-zones on a processor before swapping with z
!  rhol          : padded density after the y-lagrangian step (g cm^{-3})
!  rhoi          : padded density before the y-lagrangian step (g cm^{-3})
!  tl            : padded temperature after the y-lagrangian step (K)
!  ti            : padded temperature before the y-lagrangian step (K)
!
!    Output arguments:
!  dtime_y_hydro : minimum times step by hydro criteria
!  dt_y          : minimum hydro time step restrictions for criterion i, 
!                   of all the restrictions in angular ray (ji_ray,jk_ray) (s)
!  jdt_y         : radial zone settung minimum hydro time step restrictions
!                   for criterion i, angular zone j (s)
!
!    Include files:
!  kind_module, numerical_module
!  angular_ray_module, eos_snc_y_module, prb_cntl_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE angular_ray_module, ONLY : dyphys_c, rho_y, dt_y_state
USE eos_snc_y_module, ONLY: gam1, aesv
USE prb_cntl_module, ONLY: ihydro
USE t_cntrl_module, ONLY: tcntrl

IMPLICIT none
SAVE
!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jmin            ! inner y-array index
INTEGER, INTENT(in)              :: jmax            ! outer y-array index
INTEGER, INTENT(in)              :: ji_ray          ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray          ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: j_ray_dim       ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of radial zones on a processor before swapping with z

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(3,j_ray_dim,ik_ray_dim)            :: jdt_y         ! angular zone determining the minimum allowed time step

REAL(KIND=double), INTENT(out)                                     :: dtime_y_hydro ! minimum times step by hydro criteria
REAL(KIND=double), INTENT(out), DIMENSION(3,j_ray_dim,ik_ray_dim)  :: dt_y          ! minimum allowed time step for 

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j               ! y-array zone index

REAL(KIND=double), PARAMETER     :: dtmax = 1.d+20  ! times step without criteria
REAL(KIND=double)                :: d               ! working variable
REAL(KIND=double)                :: dmin            ! minimum value of a set of quantities
REAL(KIND=double), PARAMETER     :: frthrd = 4.d0/3.d0

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

!-----------------------------------------------------------------------
!
!        Timestep criteria.
!
!    tcntrl(4) is the fraction of a y-zone width a sonic disturbance
!     is allowed to propagate in one time step.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!            \\\\\ COURANT TIME STEP DT_Y(1,J_RAY) /////
!
!-----------------------------------------------------------------------

!----------------------------------------------------------------------!
!  Only use the Courant condition if ihydro /= 0
!----------------------------------------------------------------------!

IF ( ihydro /= 0 ) THEN

  dmin                  = dtmax

  DO j = jmin,jmax

!-----------------------------------------------------------------------
!  Find minimum sound crossing time
!-----------------------------------------------------------------------

    IF ( aesv(j,1,ji_ray,jk_ray) > zero ) THEN
      d                 = dyphys_c(j,ji_ray,jk_ray)                           &
&                       * dyphys_c(j,ji_ray,jk_ray) * rho_y(j,ji_ray,jk_ray)  &
&                       / ( DMAX1( gam1(j,ji_ray,jk_ray), frthrd ) * aesv(j,1,ji_ray,jk_ray) )
      IF ( d < dmin ) THEN
        dmin            = d
        jdt_y(1,ji_ray,jk_ray) = j
      END IF ! d < dmin
    END IF ! aesv(j,1,i_ray) > 0

  END DO

!-----------------------------------------------------------------------
!  Put common Courant condition in dt_y(1,ji_ray,jk_ray)
!-----------------------------------------------------------------------

  dt_y(1,ji_ray,jk_ray) = tcntrl(4) * DSQRT(dmin)

END IF ! ihydro /= 0

!-----------------------------------------------------------------------
!
!               \\\\\ SELECT MINIMUM TIME STEP /////
!
!    dt_y(1,j_ray) : minimum y-Currant time step
!    dt_y_state    : minimum y time step determined by state variable
!     changes
!-----------------------------------------------------------------------

dtime_y_hydro           = DMIN1( dt_y(1,ji_ray,jk_ray), dt_y_state )

RETURN
END SUBROUTINE courant_y_time_step

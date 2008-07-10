SUBROUTINE courant_z_time_step( kmin, kmax, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, dt_z, jdt_z, dtime_z_hydro )
!-----------------------------------------------------------------------
!
!    File:         courant_z_time_step
!    Module:       courant_z_time_step
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
!  kmin          : inner y-array index
!  kmax          : outer y-array index
!  ki_ray        : x (radial) index of a specific z (azimuthal) ray
!  kj_ray        : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim    : number of y (angular) zones on a processor before swapping with y
!  k_ray_dim     : number of x (radial) zones on a processor after swapping with z
!  rhol          : padded density after the y-lagrangian step (g cm^{-3})
!  rhoi          : padded density before the y-lagrangian step (g cm^{-3})
!  tl            : padded temperature after the y-lagrangian step (K)
!  ti            : padded temperature before the y-lagrangian step (K)
!
!    Output arguments:
!  dtime_z_hydro : minimum times step by hydro criteria
!  dt_z          : minimum hydro time step restrictions for criterion i, 
!                   of all the restrictions in angular ray (kj_ray,ki_ray) (s)
!  jdt_z         : radial zone settung minimum hydro time step restrictions
!                   for criterion i, angular zone k (s)
!
!    Include files:
!  kind_module, numerical_module
!  angular_ray_module, eos_snc_z_module, prb_cntl_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero

USE azimuthal_ray_module, ONLY : dzphys_c, rho_z, dt_z_state
USE eos_snc_z_module, ONLY: gam1, aesv
USE prb_cntl_module, ONLY: ihydro
USE t_cntrl_module, ONLY: tcntrl

IMPLICIT none
SAVE
!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: kmin            ! inner y-array index
INTEGER, INTENT(in)              :: kmax            ! outer y-array index
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim      ! number of y (angular) zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim       ! number of x (radial) zones on a processor after swapping with z

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(3,ij_ray_dim,k_ray_dim)            :: jdt_z         ! angular zone determining the minimum allowed time step

REAL(KIND=double), INTENT(out)                                     :: dtime_z_hydro ! minimum times step by hydro criteria
REAL(KIND=double), INTENT(out), DIMENSION(3,ij_ray_dim,k_ray_dim)  :: dt_z          ! minimum allowed time step for 

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: k               ! z-array zone index

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
!    tcntrl(7) is the fraction of a z-zone width a sonic disturbance
!     is allowed to propagate in one time step.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!            \\\\\ COURANT TIME STEP DT_Z(1,J_RAY) /////
!
!-----------------------------------------------------------------------

!----------------------------------------------------------------------!
!  Only use the Courant condition if ihydro /= 0
!----------------------------------------------------------------------!

IF ( ihydro /= 0 ) THEN

  dmin                  = dtmax

  DO k = kmin,kmax

!-----------------------------------------------------------------------
!  Find minimum sound crossing time
!-----------------------------------------------------------------------

    IF ( aesv(k,1,kj_ray,ki_ray) > zero ) THEN
      d                 = dzphys_c(k,kj_ray,ki_ray)                           &
&                       * dzphys_c(k,kj_ray,ki_ray) * rho_z(k,kj_ray,ki_ray)  &
&                       / ( DMAX1( gam1(k,kj_ray,ki_ray), frthrd ) * aesv(k,1,kj_ray,ki_ray) )
      IF ( d < dmin ) THEN
        dmin            = d
        jdt_z(1,kj_ray,ki_ray) = k
      END IF ! d < dmin
    END IF ! aesv(k,1,i_ray) > 0

  END DO

!-----------------------------------------------------------------------
!  Put common Courant condition in dt_z(1,kj_ray,ki_ray)
!-----------------------------------------------------------------------

  dt_z(1,kj_ray,ki_ray) = tcntrl(7) * DSQRT(dmin)

END IF ! ihydro /= 0

!-----------------------------------------------------------------------
!
!               \\\\\ SELECT MINIMUM TIME STEP /////
!
!-----------------------------------------------------------------------

dtime_z_hydro           = DMIN1( dt_z(1,kj_ray,ki_ray), dt_z_state )

RETURN
END SUBROUTINE courant_z_time_step

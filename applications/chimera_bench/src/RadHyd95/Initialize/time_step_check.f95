SUBROUTINE time_step_check( ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         time_step_check
!    Module:       time_step_check
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  K. R. DeNisco, Dept. of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/27/04
!
!    Purpose:
!      To compute the minimum time step as given by the Courant condition
!       and compare it with the given initial time step.
!
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ij_ray_dim    : number of y-zones on a processor before swapping
!  ik_ray_dim    : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  edit_module, evh1_global, evh1_sweep, radial_ray_module
!
!-----------------------------------------------------------------------
 
USE kind_module
USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE evh1_global, ONLY : svel, ngeomy, courant
USE evh1_sweep, ONLY : nmin, nmax, r, ei, ye, temp, gc, p

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, dtnph, &
& rho_c, t_c, ye_c, ei_c, p_c, gc_c, dx_ci, dy_ci, x_ci, ndim

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray_dim ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i          ! x-zone index
INTEGER                          :: ij_ray     ! j-index of a radial ray
INTEGER                          :: ik_ray     ! k-index of a radial ray

REAL(KIND=double)                :: width      ! x grid zone width
REAL(KIND=double)                :: widthy     ! y grid zone width
REAL(KIND=double)                :: sveln      ! reciprocal of sound crossing time
REAL(KIND=double)                :: dt_courant ! Currant time

 1001 FORMAT (' svel = 0 in subroutine time_step_check')
 1003 FORMAT (' dtnph = 0 in subroutine time_step_check')
 1005 FORMAT (' WARNING: dt_courant=',es11.3,' < dtnph=',es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

svel                   = 0.0d0
nmin                   = imin + 6
nmax                   = imax + 6

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    r   (nmin:nmax)    = rho_c(1:imax,ij_ray,ik_ray)
    ei  (nmin:nmax)    = ei_c (1:imax,ij_ray,ik_ray)
    ye  (nmin:nmax)    = ye_c (1:imax,ij_ray,ik_ray)
    temp(nmin:nmax)    = t_c  (1:imax,ij_ray,ik_ray)
    p   (nmin:nmax)    = p_c  (1:imax,ij_ray,ik_ray)
    gc  (nmin:nmax)    = gc_c (1:imax,ij_ray,ik_ray)
    IF ( ndim == 1 ) dy_ci(ij_ray) = 1.d+20
    DO i = nmin,nmax
      widthy           = dy_ci(ij_ray)
      IF ( ngeomy > 2 ) widthy = widthy * x_ci(i-6)
      width            = DMIN1( dx_ci(i-6), widthy )
      sveln            = DSQRT( gc(i) * p(i)/r(i) )/width
      svel             = DMAX1( sveln, svel )
    END DO ! i = nmin,nmax
  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

IF ( svel == zero ) THEN
  WRITE (nlog,1001)
  STOP
END IF

IF ( dtnph == zero ) THEN
  WRITE (nlog,1003)
  STOP
END IF

dt_courant             = courant/svel

IF ( dt_courant < dtnph ) THEN
  WRITE (nlog,1005) dt_courant, dtnph
END IF

RETURN
END SUBROUTINE time_step_check

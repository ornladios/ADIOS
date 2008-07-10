 SUBROUTINE set_final_angular_grid
!-----------------------------------------------------------------------
!
!    File:         set_final_angular_grid
!    Module:       set_final_angular_grid
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/05/05
!
!    Purpose:
!      To Set the angular grid at the end of the hydro step.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  i_ray_dim : number of radial rays on a processor
!  nx        : x-array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  cycle_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half, epsilon

USE cycle_module, ONLY : ncycle
USE radial_ray_module, ONLY : jmin, jmax, y_ei, dy_ci, y_ci, y_ef, dy_cf, &
& y_cf, y_shft, ncy_shift, dy_shift, tb_dy_shift, time, t_bounce

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                               :: first_y = .true.

REAL(KIND=double)                     :: y_shift          ! y-grid shift parameter

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!    \\\\\ Y-GRID SHIFT PARAMETER FOR ADDITIONAL DISSIPATION /////
!
!-----------------------------------------------------------------------

IF ( first_y ) THEN
  first_y                 = .false.
  y_shift                 = zero
  IF ( jmax - jmin == 1 ) THEN
    y_shift               = dy_shift * half * ( y_ei(jmax+1) - y_ei(jmin) )
  END IF ! jmax - jmin = 1
  IF ( jmax - jmin > 1 ) THEN
    y_shift               = dy_shift * dy_ci(jmax/2)
  END IF ! jmax - jmin > 1
END IF ! first_y

!-----------------------------------------------------------------------
!
!                         \\\\\ EULERIAN /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set final grid to the initial grid
!-----------------------------------------------------------------------

y_ef                      = y_ei
dy_cf                     = dy_ci
y_cf                      = y_ci

!-----------------------------------------------------------------------
!  If y_shift == 'ye', wiggle grid for additional lateral dissipation;
!   recalculate solid angles
!-----------------------------------------------------------------------

IF ( y_shft == 'ye' ) THEN
  IF ( ncycle > ncy_shift  .and.  jmax - jmin > 0 ) THEN
    y_ef(jmin+1:jmax)     = y_ei(jmin+1:jmax) + ( 2 * MOD( ncycle, 2 ) - 1 ) * y_shift
    y_cf(jmin:jmax)       = half * ( y_ef(jmin:jmax) + y_ef(jmin+1:jmax+1) )
    dy_cf(jmin)           = y_ef(jmin+1) - y_ef(jmin)
    dy_cf(jmax)           = y_ef(jmax+1) - y_ef(jmax)
    IF ( t_bounce > epsilon             .and.  &
&        time - t_bounce > tb_dy_shift  .and.  &
&        MOD( ncycle, 2 ) == 0 ) y_shft = 'no'
    END IF ! ncycle > ncy_shift  .and.  jmax - jmin > 0
  CALL solid_angle
END IF !  y_shft == 'ye'

RETURN
END SUBROUTINE set_final_angular_grid

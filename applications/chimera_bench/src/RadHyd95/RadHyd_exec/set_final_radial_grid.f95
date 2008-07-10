 SUBROUTINE set_final_radial_grid( ij_ray_dim, ik_ray_dim, nx )
!-----------------------------------------------------------------------
!
!    File:         set_final_radial_grid
!    Module:       set_final_radial_grid
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/05/05
!
!    Purpose:
!      To set the radial grid at the end of the hydro step.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ij_ray_dim : number of y-zones on a processor before swapping (used in MPI version)
!  ik_ray_dim : number of z-zones on a processor before swapping (used in MPI version)
!  nx         : x-array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  cycle_module, edit_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half, frpi, epsilon

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nprint, nlog
USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax,  m_grid,   &
& u_e, x_ei, x_ei, dx_ci, x_ci, x_ef, dx_cf, x_cf, x_el, dx_cl, x_cl, dtnph, &
& d_omega, lagr

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER                               :: ij_ray_dim    ! number of y-zones on a processor before swapping (used in MPI version)
INTEGER                               :: ik_ray_dim    ! number of z-zones on a processor before swapping (used in MPI version)
INTEGER                               :: nx            ! x-array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                               :: first_y = .true.

INTEGER                               :: i             ! x-array index

REAL(KIND=double)                     :: y_shift       ! y-grid shift parameter

REAL(KIND=double), DIMENSION(nx)      :: u_grid        ! x-velocity of the grid

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                   \\\\\ LAGRANGIAN GRID /////
!
!-----------------------------------------------------------------------

IF ( lagr == 'ye' ) THEN

  x_ef (imin:imax+1)      = x_el (imin:imax+1,1,1)
  dx_cf(imin:imax)        = dx_cl(imin:imax  ,1,1)
  x_cf (imin:imax)        = x_cl (imin:imax  ,1,1)

  RETURN

END IF ! lagr == 'ye'

!-----------------------------------------------------------------------
!
!                 \\\\\ EULERIAN, MOVING GRID OFF /////
!
!-----------------------------------------------------------------------

IF ( lagr == 'no'  .and.  m_grid == 'no' ) THEN

!-----------------------------------------------------------------------
!  Set final grid to the initial grid
!-----------------------------------------------------------------------

  x_ef                    = x_ei
  dx_cf                   = dx_ci
  x_cf                    = x_ci

  RETURN

END IF ! lagr == 'no'  .and.  m_grid == 'no'

!-----------------------------------------------------------------------
!
!                      \\\\\ MOVING GRID ON /////
!
!-----------------------------------------------------------------------

IF ( lagr == 'no'  .and.  m_grid == 'ye' ) THEN

!-----------------------------------------------------------------------
!  Compute final grip position
!-----------------------------------------------------------------------

  DO i = imin,imax+1
    u_grid(i)            = SUM( u_e(i,jmin:jmax,kmin:kmax) * d_omega(jmin:jmax,kmin:kmax) )/frpi
  END DO ! i

  x_ef (imin:imax+1)     = x_ei(imin:imax+1) + u_grid(imin:imax+1) * dtnph
  dx_cf(imin:imax)       = x_ef(imin+1:imax+1) - x_ef(imin:imax)
  x_cf (imin:imax)       = half * ( x_ef(imin+1:imax+1) + x_ef(imin:imax) )

  RETURN

END IF ! lagr == 'no'  .and.  m_grid == 'ye'

RETURN
END SUBROUTINE set_final_radial_grid

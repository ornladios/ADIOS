SUBROUTINE editMDgy( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_in, &
& y_in, z_in, time, t_tb, grav_y_cMD, nprint )
!-----------------------------------------------------------------------
!
!    File:         editMDgy
!    Module:       editMDgy
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/08/07
!
!    Purpose:
!  To perform an edit of the x-gravitation al force
!
!    Subprograms call:
!  date_and_time_print : fetches and records the date and the time
!
!    Input arguments:
!  imin         : minimum x-array index for the edit
!  imax         : maximum x-array index for the edit
!  nx           : x-array extent
!  jmin         : minimum y-array index for the edit
!  jmax         : maximum y-array index for the edit
!  ny           : y_array extent
!  x_in         : radial midpoint of zone (cm)
!  y_in         : angular midpoint of zone
!  time         : elapsed time (s)
!  t_tb         : time from core bounce (s)
!  grav_y_cMD   : zone-centered y-component of gravitational acceleration (cm s^{-2} g^{-1})
!  iunitp       : unit number to print models
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                :: imin           ! minimum x-array index for the edit
INTEGER, INTENT(in)                                :: imax           ! maximum x-array index for the edit
INTEGER, INTENT(in)                                :: nx             ! x-array extent

INTEGER, INTENT(in)                                :: jmin           ! minimum y-array index for the edit
INTEGER, INTENT(in)                                :: jmax           ! maximum y-array index for the edit
INTEGER, INTENT(in)                                :: ny             ! y_array extent

INTEGER, INTENT(in)                                :: kmin           ! minimum z-array index for the edit
INTEGER, INTENT(in)                                :: kmax           ! maximum z-array index for the edit
INTEGER, INTENT(in)                                :: nz             ! z_array extent

INTEGER, INTENT(in)                                :: nprint         ! unit number to print models

REAL(KIND=double), INTENT(in)                      :: time           ! elapsed time (s)
REAL(KIND=double), INTENT(in)                      :: t_tb           ! time from core bounce (s)
REAL(KIND=double), INTENT(in), DIMENSION(nx)       :: x_in           ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny)       :: y_in           ! angular midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nz)       :: z_in           ! z (azimuthal) midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,ny,nz) :: grav_y_cMD     ! zone-centered y-component of gravitational acceleration [cm s^{-2} g^{-1}]

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

WRITE (nprint) imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz
WRITE (nprint) time, t_tb
WRITE (nprint) x_in
WRITE (nprint) y_in
WRITE (nprint) z_in
WRITE (nprint) grav_y_cMD

RETURN
END SUBROUTINE editMDgy
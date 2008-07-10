SUBROUTINE editMDnl( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, n, &
& nnu, x_in, y_in, z_in, time, t_tb, lumMD, nprint )
!-----------------------------------------------------------------------
!
!    File:         editMDnl
!    Module:       editMDnl
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/22/05
!
!    Purpose:
!  To perform a 2-D edit of the neutrino luminosities.
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
!  kmin         : minimum z-array index for the edit
!  kmax         : maximum z-array index for the edit
!  nz           : z dimension
!  n            : neutrino flavor index
!  nnu          : neutrino flavor array extent
!  x_in         : radial midpoint of zone (cm)
!  y_in         : angular midpoint of zone
!  z_in         : z (azimuthal) midpoint of zone
!  time         : elapsed time (s)
!  t_tb         : time from core bounce (s)
!  lumMD        : neutrino luminosity [B]
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

INTEGER, INTENT(in)                                    :: imin        ! minimum x-array index for the edit
INTEGER, INTENT(in)                                    :: imax        ! maximum x-array index for the edit
INTEGER, INTENT(in)                                    :: nx          ! x-array extent

INTEGER, INTENT(in)                                    :: jmin        ! minimum y-array index for the edit
INTEGER, INTENT(in)                                    :: jmax        ! maximum y-array index for the edit
INTEGER, INTENT(in)                                    :: ny          ! y_array extent

INTEGER, INTENT(in)                                    :: kmin           ! minimum z-array index for the edit
INTEGER, INTENT(in)                                    :: kmax           ! maximum z-array index for the edit
INTEGER, INTENT(in)                                    :: nz             ! z_array extent

INTEGER, INTENT(in)                                    :: n           ! neutrino flavor index
INTEGER, INTENT(in)                                    :: nnu         ! neutrino flavor array extent

INTEGER, INTENT(in)                                    :: nprint      ! unit number to print models

REAL(KIND=double), INTENT(in)                          :: time        ! elapsed time (s)
REAL(KIND=double), INTENT(in)                          :: t_tb        ! time from core bounce (s)
REAL(KIND=double), INTENT(in), DIMENSION(nx)           :: x_in        ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny)           :: y_in        ! angular midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nz)           :: z_in           ! z (azimuthal) midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,ny,nz) :: lumMD       ! neutrino luminosity [B]

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

WRITE (nprint) imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz
WRITE (nprint) time, t_tb
WRITE (nprint) x_in
WRITE (nprint) y_in
WRITE (nprint) z_in
WRITE (nprint) lumMD(:,n,:,:)

RETURN
END ! SUBROUTINE editMDnl

MODULE conditions
!===============================================================================
!  This module contains data on the current time and thermodynamic conditions.
!===============================================================================

USE kind_module, ONLY : double

REAL(KIND=double) :: t        ! Time at the beginning of the current timestep
REAL(KIND=double) :: tt       ! Trial time for end of current timestep
REAL(KIND=double) :: tdel     ! Trial duration of timestep
REAL(KIND=double) :: t9t      ! Temperature (GK) at trial time
REAL(KIND=double) :: rhot     ! Density (g/cc) at trial time

END MODULE conditions

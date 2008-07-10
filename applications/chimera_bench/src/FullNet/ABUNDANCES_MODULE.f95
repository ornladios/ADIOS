MODULE abundances
!===============================================================================
!  This module contains the abundances of the nuclear species, at the previous 
!  time (yo), current time (y), and trial time (yt), as well as the time 
!  derivatives (ydot), at the trial time, and the abundance change due to the 
!  Newton-Raphson iteration.  The size of these arrays is allocated in the 
!  external routine, where the initial values are set.
!===============================================================================

USE kind_module, ONLY : double
USE nuc_number

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: yo
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: yt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ydot
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dy

END MODULE abundances

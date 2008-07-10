!______________________________________________________________________
! module for TOV potential
!
!______________________________________________________________________

MODULE tov_potential_module

USE kind_module
SAVE

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: effpot_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: effpot_e

END MODULE tov_potential_module

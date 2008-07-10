MODULE part_funct_data
!-----------------------------------------------------------------------------
!  This module contains the nuclear partition function data.  The value of
!  the partition function is interpolated over a grid of temperatures, t9i.
!  g is the interpolation data (24,ny), gg is the current parition function,
!  and angm is the J. gg(0) and angm(0) are placeholders for non-nuclei
!  The array sizes are set in read_nuclear_data.
!-----------------------------------------------------------------------------

USE kind_module, ONLY : double

REAL(KIND=double)                              :: t9i(24)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: gg,angm
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: g

END MODULE part_funct_data

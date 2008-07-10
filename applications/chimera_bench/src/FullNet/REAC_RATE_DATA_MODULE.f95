MODULE reac_rate_data
!-----------------------------------------------------------------------------
!  This module contains the data necessary to calculate the reaction rates
!  and to map to each species those reactions which affect it.
!-----------------------------------------------------------------------------

USE kind_module, ONLY : double

INTEGER                                      :: nan(3)
INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: la,le
INTEGER, ALLOCATABLE, DIMENSION(:)           :: mu1,mu2,mu3
INTEGER, ALLOCATABLE, DIMENSION(:)           :: n11,n21,n22,n31,n32,n33

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: a1,a2,a3,b1,b2,b3 ! dim(nan(i))

END MODULE reac_rate_data


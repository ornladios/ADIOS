MODULE ffn_data
!-----------------------------------------------------------------------------
!  This module contains the data to calculate weak reactions according to
!  Fuller, Fowler, Neuman (1982,1985).
!-----------------------------------------------------------------------------

USE kind_module, ONLY : double

INTEGER                                        :: nffn           ! The number of FFN reactions

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rf,r1,r2       ! dim(nffn)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: ffnsum,ffnenu  ! dim(nffn,143)

END MODULE ffn_data

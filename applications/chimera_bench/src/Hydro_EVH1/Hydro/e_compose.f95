SUBROUTINE e_compose( nmin, nmax )

!-----------------------------------------------------------------------
!  Calculates the total energy per gram (not including gravity), e(n), 
!   from the internal and kinetic energies
!-----------------------------------------------------------------------

USE evh1_sweep, ONLY: e, ei, ekin
     
IMPLICIT none

INTEGER, INTENT(in)                    :: nmin     ! minimum padded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum padded array index

!-----------------------------------------------------------------------

CALL e_kinetic( nmin, nmax )

e(nmin:nmax) = ei(nmin:nmax) + ekin(nmin:nmax)

RETURN
END SUBROUTINE e_compose

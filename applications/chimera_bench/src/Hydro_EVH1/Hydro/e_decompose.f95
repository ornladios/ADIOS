SUBROUTINE  e_decompose( nmin, nmax )
!-----------------------------------------------------------------------
!  Calculates the internal energy per gram, ei(n), from the total energy
!   (not including gravity) and the kinetic energies
!-----------------------------------------------------------------------


USE evh1_sweep, ONLY: e, ei, ekin
     
IMPLICIT none

INTEGER, INTENT(in)                    :: nmin     ! minimum padded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum padded array index

!-----------------------------------------------------------------------

CALL e_kinetic( nmin, nmax )

ei(nmin:nmax) = e(nmin:nmax) - ekin(nmin:nmax)

RETURN
END SUBROUTINE e_decompose

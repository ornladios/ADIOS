SUBROUTINE  e_decompose_v( nmin, nmax )
!-----------------------------------------------------------------------
!  Calculates the partial energy per gram, e_v(n), from the total and
!  kinetic energies
!-----------------------------------------------------------------------


USE evh1_sweep, ONLY: e, e_v, ekin
     
IMPLICIT none

INTEGER, INTENT(in)                    :: nmin     ! minimum padded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum padded array index

!-----------------------------------------------------------------------

CALL e_kinetic( nmin, nmax )

e_v(nmin:nmax) = e(nmin:nmax) - ekin(nmin:nmax)

RETURN
END SUBROUTINE e_decompose_v

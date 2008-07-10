SUBROUTINE e_compose_v( nmin, nmax )

!-----------------------------------------------------------------------
!  Calculates the total energy per gram (including gravity), e(n),  from
!   the partial energy and kinetic energies
!-----------------------------------------------------------------------

USE evh1_sweep, ONLY: e, e_v, ekin
     
IMPLICIT none

INTEGER, INTENT(in)                    :: nmin     ! minimum padded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum padded array index

!-----------------------------------------------------------------------

CALL e_kinetic( nmin, nmax )

e(nmin:nmax) = e_v(nmin:nmax) + ekin(nmin:nmax)

RETURN
END SUBROUTINE e_compose_v

SUBROUTINE e_compose_g( nmin, nmax )

!-----------------------------------------------------------------------
!  To calculate the total energy per gram (including gravity), e(n),
!   from the internal, kinetic, and gravitational potential energies
!-----------------------------------------------------------------------

USE evh1_zone, ONLY: imax
USE evh1_sweep, ONLY: e, ei, egrav, ekin
     
IMPLICIT none

INTEGER, INTENT(in)                    :: nmin     ! minimum padded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum padded array index

!-----------------------------------------------------------------------

CALL e_kinetic( nmin, nmax )

e(nmin:nmax) = ei(nmin:nmax) + ekin(nmin:nmax) + egrav(nmin:nmax)

RETURN
END SUBROUTINE e_compose_g

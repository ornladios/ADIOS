SUBROUTINE  e_decompose_g( nmin, nmax )
!-----------------------------------------------------------------------
!  To calculate the internal energy per gram, ei(n), from the total 
!   energy (including gravit) by subrtracting off the kinetic and 
!   gravitational potential energies
!-----------------------------------------------------------------------

USE evh1_global, ONLY: small
USE evh1_sweep, ONLY: e, ei, egrav, ekin, degrav
     
IMPLICIT none

INTEGER, INTENT(in)                    :: nmin     ! minimum padded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum padded array index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CALL e_kinetic( nmin, nmax )

ei(nmin:nmax) = e(nmin:nmax) - ekin(nmin:nmax) - egrav(nmin:nmax) + degrav(nmin:nmax)

RETURN
END SUBROUTINE e_decompose_g

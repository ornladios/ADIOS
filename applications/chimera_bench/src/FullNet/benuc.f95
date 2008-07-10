SUBROUTINE benuc( y, enb, enm, ytot, ztot, atot )
!===============================================================================
!  This routine finds moments of the abundance distribution useful for
!  hydrodynamics, including the total abundance, electron fraction, binding 
!  energy, and mass excess energy and outputs in and mol/g and ergs/g. 
!===============================================================================

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one

USE constants, ONLY : epmev, avn
USE nuclear_data, ONLY : zz, aa, be
USE nuc_number, ONLY : ny

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: y(ny)  ! abundance fractions

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out) :: enb      ! mean binding energy per particle
REAL(KIND=double), INTENT(out) :: enm      ! rest mass energy
REAL(KIND=double), INTENT(out) :: ytot     ! abundance fraction
REAL(KIND=double), INTENT(out) :: ztot     ! electron fraction
REAL(KIND=double), INTENT(out) :: atot     ! mass fraction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER   :: mex_p=7.28899d0, mex_n=8.07144d0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ytot         = SUM(y)
ztot         = SUM(zz*y)
atot         = SUM(aa*y)
enb          = SUM(be*y)

enm          = mex_p * ztot + mex_n * ( atot - ztot ) - enb

!-----------------------------------------------------------------------
!  Change units from MeV/nucleon to erg/g
!-----------------------------------------------------------------------

enb          = epmev * avn * enb
enm          = epmev * avn * enm

RETURN
END SUBROUTINE benuc

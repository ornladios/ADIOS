SUBROUTINE eos_y_reset( jmin, jmax, rho, t, ye, ji_ray, jk_ray, &
& reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         eos_y_reset
!    Module:       eos_y_reset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/17/05
!
!    Purpose:
!      To reset the thermodymamic and neutrino rate tables.
!
!    Subprograms called:
!  esrgnz_comp_y  : recomputes (always) EOS values at cube edges
!  esrgnz_y       : recomputes (if necessary) EOS values at cube edges
!  eqstz_y        : interpolates EOS values
!  gammaz_y       : computes EOS gammas
!
!    Input arguments:
!  jmin           : minimum y-array index
!  jmax           : maximum y-array index
!  rho            : angular array of matter density (g cm{-3})
!  t              : angular array of matter temperature (K)
!  ye             : angular array of matter electron fraction
!  ji_ray         : x (radial) index of a specific angular ray
!  jk_ray         : z (azimuthal) index of a specific angular ray
!  reset_comp_eos : composition EOS reset flag
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : ny

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos  ! composition EOS reset flag

INTEGER, INTENT(in)              :: jmin            ! minimum y-array index
INTEGER, INTENT(in)              :: jmax            ! maximum y-array index
INTEGER, INTENT(in)              :: ji_ray          ! x (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray          ! z (azimuthal) index of a specific angular ray

REAL(KIND=double), INTENT(in), DIMENSION(ny) :: rho ! angular array of matter density array (g cm{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: t   ! angular array of matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: ye  ! angular array of matter matter electron fraction array

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ RESET EOS TABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Recompute thermodynamic quantities
!-----------------------------------------------------------------------

IF ( reset_comp_eos == 'ye' ) THEN
  CALL esrgnz_comp_y( jmin, jmax, rho, t, ye, ji_ray, jk_ray )
ELSE
  CALL esrgnz_y( jmin, jmax, rho, t, ye, ji_ray, jk_ray )
END IF
CALL eqstz_y( jmin, jmax, rho, t, ye, ji_ray, jk_ray )
CALL gammaz_y( jmin, jmax, rho, t, ji_ray, jk_ray )

RETURN
END SUBROUTINE eos_y_reset

SUBROUTINE eos_z_reset( kmin, kmax, rho, t, ye, ki_ray, kj_ray, &
& reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         eos_z_reset
!    Module:       eos_z_reset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To reset the thermodymamic and neutrino rate tables.
!
!    Subprograms called:
!  esrgnz_comp_z  : recomputes (always) EOS values at cube edges
!  esrgnz_z       : recomputes (if necessary) EOS values at cube edges
!  eqstz_z        : interpolates EOS values
!  gammaz_z       : computes EOS gammas
!
!    Input arguments:
!  kmin           : minimum z-array index
!  kmax           : maximum z-array index
!  rho            : azimuthal array of matter density (g/cm**3)
!  t              : azimuthal array of matter temperature (K)
!  ye             : azimuthal array of matter electron fraction
!  ki_ray         : x (radial) index of a specific z (azimuthaal) ray
!  kj_ray         : y (angular) index of a specific z (azimuthaal) ray
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
USE array_module, ONLY : nz

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos  ! composition EOS reset flag

INTEGER, INTENT(in)              :: kmin            ! minimum z-array index
INTEGER, INTENT(in)              :: kmax            ! maximum z-array index
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthaal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthaal) ray

REAL(KIND=double), INTENT(in), DIMENSION(nz) :: rho ! azimuthal array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: t   ! azimuthal array of matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: ye  ! azimuthal array of matter matter electron fraction array

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
  CALL esrgnz_comp_z( kmin, kmax, rho, t, ye, ki_ray, kj_ray )
ELSE
  CALL esrgnz_z( kmin, kmax, rho, t, ye, ki_ray, kj_ray )
END IF
CALL eqstz_z( kmin, kmax, rho, t, ye, ki_ray, kj_ray )
CALL gammaz_z( kmin, kmax, rho, t, ki_ray, kj_ray )

RETURN
END SUBROUTINE eos_z_reset

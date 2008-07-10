SUBROUTINE eos_z( ki_ray, kj_ray )

!***********************************************************************

USE kind_module, ONLY : double
USE eos_bck_module, ONLY : jshel, dbck, tbck, yebck, erad, prad, enu, pnu, sneu, &
& dtran, etot, eh, ee, ed, egy0, ptotbck, ph, pe, pd, stotbck, sh, se, sd, dedt
USE eos_snc_z_module, ONLY : nse
USE physcnst_module, ONLY : dmnp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (angular) index of a specific z (azimuthal) ray

REAL(KIND=double), PARAMETER     :: asig = 8.56d-8
REAL(KIND=double)                :: at4           ! photon energy (MeV fm^{-3})
REAL(KIND=double)                :: srad          ! entropy (nucleon^{-1})
REAL(KIND=double)                :: yebck_save    ! saved value of yebck
REAL(KIND=double)                :: dtest         ! estimate of transition density to nuclear matter

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!*************************************************** radiation *********

at4             = asig * tbck * tbck * tbck * tbck
erad            = at4/dbck
prad            = at4/3.d0
srad            = ( erad + prad/dbck )/tbck
enu             = erad
pnu             = prad
sneu            = srad

!************************************************** electrons *********
!*************************************************no neutrinos ********

CALL lectron

!**************************************************** nucleons *********

CALL eos0

!........Store yebck

yebck_save      = yebck

IF ( nse(jshel,kj_ray,ki_ray) == 1 ) THEN

!........Matter is in NSE

  yebck         = DMIN1( yebck, 0.49999d0 )

  dtest         = .16d0 * ( 1.d0 - 3.d0 * ( .5d0 - yebck )**2 )
  IF ( dtran <= 0.0d0 ) dtran = dtest
  dtran         = DMIN1( dtran, dtest )

!........Transition to nuclear matter

  IF ( dbck >= dtran ) THEN
    CALL eosnm
  ELSE
    CALL saha
  END IF

!........Matter is not in NSE

ELSE
  CALL eosnuc_z
END IF

!************************************************ get totals ***********
!-----------------------------------------------------------------------
!        Energy offset
!
!  The subroutine  eosnm, saha, and eosnuc_x return the energy including
!   the binding energy.
!  The quantity - yebck * dmnp is added so that the binding energy is
!   relative to free neutrons rather than free neutrons and protons
!  An energy offset 8.9d0 MeV/nucleon is added to prevent the internal
!   energy from becoming negative
!-----------------------------------------------------------------------

egy0            = -8.9d0
etot            = eh + erad + ee + ed - egy0 - yebck * dmnp
ptotbck         = ph + prad + pe + pd
stotbck         = sh + srad + se + sd

!***************************************** approx for energy inversion *

dedt            = se + sh + erad*4.d0/tbck + ed/tbck
yebck           = yebck_save

RETURN
END SUBROUTINE eos_z

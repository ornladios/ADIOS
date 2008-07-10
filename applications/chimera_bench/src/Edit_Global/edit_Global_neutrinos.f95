SUBROUTINE edit_Global_neutrinos( imin, imax, nx, nnu, time, t_tb, x_c,  &
&  nse_min, nse_max, lum_mean, lum_min, lum_max, sigma_lum, e_rms_mean, &
&  e_rms_min, e_rms_max, sigma_e_rms, c_shock, nprint )
!-----------------------------------------------------------------------
!
!    File:         edit_Global_neutrinos
!    Module:       edit_Global_neutrinos
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/01/07
!
!    Purpose:
!      To edit the minimum, maximum, and mean values of neutrino luminosities
!       amd rms energies as a function of the radius.
!
!    Subprograms called:
!  date_and_time_print : fetches and prints the date and time
!
!    Input arguments:
!  imin         : minimum x-array index for the edit
!  imax         : maximum x-array index for the edit
!  nx           : x-array extent
!  nnu          : neutrino flavor array extent
!  time         : elapsed time
!  t_tb         : time from core bounce
!  x_c          : radial midpoint of zone (cm)
!  nse_min      : minimum value of NSE-nonNSE boundary
!  nse_max      : maximum value of NSE-nonNSE boundary
!  lum_mean     : mass averaged n-neutrino luminosity (foes s^{1})
!  lum_min      : minimum n-neutrino luminosity (foes s^{1}
!  lum_max      : maximum  n-neutrino luminosity (foes s^{1})
!  sigma_lum    : RMS neutrino n-neutrino luminosity (foes s^{1})
!  e_rms_mean   : mass averaged neutrino n-neutrino rms energy (MeV)
!  e_rms_min    : minimum neutrino n-neutrino rms energy (MeV)
!  e_rms_max    : maximum neutrino n-neutrino rms energy (MeV)
!  sigma_e_rms  : RMS neutrino n-neutrino rms energy (MeV))
!  c_shock      : shock location
!  nprint       : unit numberto print data
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module
!  edit_head
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE edit_module, ONLY : head

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                              :: imin           ! minimum unshifted radial zone index
INTEGER, INTENT(in)                              :: imax           ! maximum unshifted radial zone index
INTEGER, INTENT(in)                              :: nx             ! x-array extent
INTEGER, INTENT(in)                              :: nnu            ! neutrino flavor array extent
INTEGER, INTENT(in)                              :: nse_min        ! minimum value of NSE-nonNSE boundary
INTEGER, INTENT(in)                              :: nse_max        ! maximum value of NSE-nonNSE boundary
INTEGER, INTENT(in)                              :: nprint         ! unit numberto print data

CHARACTER (len=2), INTENT(in), DIMENSION(nx)     :: c_shock        ! shock location

REAL(KIND=double), INTENT(in)                    :: time           ! elapsed time
REAL(KIND=double), INTENT(in)                    :: t_tb           ! time from core bounce
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: x_c            ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: lum_mean       ! mass averaged n-neutrino luminosity (foes s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: lum_min        ! minimum n-neutrino luminosity (foes s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: lum_max        ! maximum  n-neutrino luminosity (foes s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: sigma_lum      ! RMS neutrino n-neutrino luminosity (foes s^{1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: e_rms_mean     ! mass averaged neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: e_rms_min      ! minimum neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: e_rms_max      ! maximum neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu) :: sigma_e_rms    ! RMS neutrino n-neutrino rms energy (MeV))

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=2), DIMENSION(nx)                 :: c_nse

INTEGER                                          :: i             ! do index

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1a)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (9x,'Elapsed time=',es14.7,10x,' Time from bounce=',es14.7/)
    9 FORMAT (41x,'Neutrino Luminosities and RMS energies')
   11 FORMAT (40x,40('-')/)
  101 FORMAT ('   j          r       enu <lum> enu lum mn enu lum mx sigma_elum &
& anu <lum> anu lum mn anu lum mx sigma_alum  xnu <lum> xnu lum mn xnu lum mx&
& sigma_xlum  znu <lum> znu lum mn znu lum mx sigma_zlum  enu <rms> enu rms mn&
& enu rms mx sigma_erms  anu <rms> anu rms mn anu rms mx sigma_arms  xnu <rms>&
& xnu rms mn xnu rms mx sigma_xrms  znu <rms> znu rms mn znu rms mx sigma_zrms'/)
  103 FORMAT (1x,i4,1x,a2,es11.3,a2,32es11.3)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                          \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

DO i = imin,imax+1
  IF ( i < nse_min  .and.  i < nse_max ) THEN
    c_nse(i)            = '  '
  ELSE IF ( i < nse_max ) THEN
    c_nse(i)            = '* '
  ELSE
    c_nse(i)            = '**'
  END IF
END DO ! i

!-----------------------------------------------------------------------
!  Neutrino luminosities and rms energies
!-----------------------------------------------------------------------

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( nprint )
WRITE (nprint,7) time,t_tb
WRITE (nprint,9)
WRITE (nprint,11)
WRITE (nprint,101)

DO i = imax+1, imin, -1
  WRITE (nprint,103) i,c_nse(i),x_c(i),c_shock(i), lum_mean(i,1), lum_min(i,1), &
&  lum_max(i,1), sigma_lum(i,1), lum_mean(i,2), lum_min(i,2), lum_max(i,2), &
&  sigma_lum(i,2), lum_mean(i,3), lum_min(i,3), lum_max(i,3), sigma_lum(i,3), &
&  lum_mean(i,4), lum_min(i,4), lum_max(i,4), sigma_lum(i,4), e_rms_mean(i,1), &
&  e_rms_min(i,1), e_rms_max(i,1), sigma_e_rms(i,1), e_rms_mean(i,2), &
&  e_rms_min(i,2), e_rms_max(i,2), sigma_e_rms(i,2), e_rms_mean(i,3), &
&  e_rms_min(i,3), e_rms_max(i,3), sigma_e_rms(i,3), e_rms_mean(i,4), &
&  e_rms_min(i,4), e_rms_max(i,4), sigma_e_rms(i,4)
END DO ! i

RETURN
END SUBROUTINE edit_Global_neutrinos

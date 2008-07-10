SUBROUTINE eqstt_y( i, j, ji_ray, jk_ray, rho, t, ye, esv, esvd, esvt, &
& esvy )
!-----------------------------------------------------------------------
!
!    File:         eqstt_y
!    Module:       eqstt_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/18/05
!
!    Purpose:
!      To interpolate the logs of equation of state quantities from a
!       local table of nearest entries created for each radial shell
!       midpoint, and to compute their derivatives with respect to
!       rho, t, and ye from the interpolation formulae.
!
!    Variables that must be passed through common:
!
!  estble(i,j,ida,ita,iya,ji_ray,jk_ray) : equation of state table array.
!
!    Subprograms called:
!  eosidx_y : computes EOS table indices and interpolation coefficients
!
!    Input arguments:
!  i     : thermodynamic function index.
!   i = 1 : pressure.
!   i = 2 : energy.
!   i = 3 : entropy.
!   i = 4 : neutron chemical potential.
!   i = 5 : proton chemical potential.
!   i = 6 : electron chemical potential.
!   i = 7 : free neutron mass fraction.
!   i = 8 : free proton mass fraction.
!   i = 9 : heavy nucleus mass fraction.
!   i = 10: heavy nucleus mass number.
!   i = 11: heavy nucleus charge number.
!   i = 12: gamma1.
!  j        : y (angular) zone index.
!  ji_ray   : x (radial) index of a specific angular ray
!  jk_ray   : z (azimuthal) index of a specific angular ray
!  rho      : matter density (g/cm**3).
!  t        : matter temperature (K).
!  ye       : matter electron fraction.
!
!    Output arguments:
!  esv      : interpolated value of thermodynamic variable i.
!  esvd     : derivative of esv wrt rho.
!  esvt     : derivative of esv wrt t.
!  esvy     : derivative of esv wrt ye.
!
!    Include files:
!  kind_module, numerical_module,
!  edit_module, eos_snc_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon

USE edit_module, ONLY : nprint, nlog
USE eos_snc_y_module, ONLY : escnst, estble

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: i             ! eos parameter index
INTEGER, INTENT(in)              :: j             ! y (angular) zone index.
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific angular ray

REAL(KIND=double), INTENT(in)    :: rho           ! density (g cm^-3)
REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: esv           ! interpolated eos parameter
REAL(KIND=double), INTENT(out)   :: esvd          ! d(esv)/d(rho)
REAL(KIND=double), INTENT(out)   :: esvt          ! d(esv)/d(t)
REAL(KIND=double), INTENT(out)   :: esvy          ! d(esv)/d(ye)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER, PARAMETER               :: ida = 1       ! lower cube table density index
INTEGER, PARAMETER               :: idap1 = 2     ! upper cube table density index
INTEGER, PARAMETER               :: ita = 1       ! lower cube table temperature index
INTEGER, PARAMETER               :: itap1 = 2     ! upper cube table temperature index
INTEGER, PARAMETER               :: iya = 1       ! lower cube table electron fraction index
INTEGER, PARAMETER               :: iyap1 = 2     ! upper cube table electron fraction index

REAL(KIND=double)                :: fdp           ! position of rho in grid wrt lower cube index
REAL(KIND=double)                :: fdm           ! position of rho in grid wrt upper cube index
REAL(KIND=double)                :: fdd           ! d(fd)/d(rho)
REAL(KIND=double)                :: ftp           ! position of t in grid wrt lower cube index
REAL(KIND=double)                :: ftm           ! position of t in grid wrt upper cube index
REAL(KIND=double)                :: ftt           ! d(ft)/d(t)
REAL(KIND=double)                :: fyp           ! position of ye in grid wrt lower cube index
REAL(KIND=double)                :: fym           ! position of ye in grid wrt upper cube index
REAL(KIND=double)                :: fyy           ! d(fy)/d(ye)

REAL(KIND=double)                :: ln_10         ! ln(10)
REAL(KIND=double), PARAMETER     :: alogmx = 2.d2 ! warn if |log10(esv)| > alogmx
REAL(KIND=double), PARAMETER     :: realmx = 1.d200 ! warn if |esv| > realmxx

REAL(KIND=double)                :: av111         ! scalar table entry for interpolation
REAL(KIND=double)                :: av211         ! scalar table entry for interpolation
REAL(KIND=double)                :: av121         ! scalar table entry for interpolation
REAL(KIND=double)                :: av112         ! scalar table entry for interpolation
REAL(KIND=double)                :: av221         ! scalar table entry for interpolation
REAL(KIND=double)                :: av212         ! scalar table entry for interpolation
REAL(KIND=double)                :: av122         ! scalar table entry for interpolation
REAL(KIND=double)                :: av222         ! scalar table entry for interpolation

REAL(KIND=double)                :: besv          ! log10(esvtmp)
REAL(KIND=double)                :: besvd         ! d(besv)/d(rho)
REAL(KIND=double)                :: besvt         ! d(besv)/d(t)
REAL(KIND=double)                :: besvy         ! d(besv)/d(ye)
REAL(KIND=double)                :: esvtmp        ! esv + escnst(i,j,ji_ray,jk_ray) + epsilon

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (1x)
 1003 FORMAT (1x,20x,'***parameter out of range in eqstt_y***')
 1005 FORMAT (' j=',i4,' ji_ray=',i4,' jk_ray=',i4,' rho=',1pe10.3,' t=',1pe10.3, &
& ' ye=',1pe10.3)
 1007 FORMAT (' i=',i4,' esv=',1pe10.3)
 1009 FORMAT (' i=',i4,' esvd=',1pe10.3)
 1011 FORMAT (' i=',i4,' esvt=',1pe10.3)
 1013 FORMAT (' i=',i4,' esvy=',1pe10.3)
 1015 FORMAT (' av111,...=',8(1pe11.3))
 1017 FORMAT (' fdm,ftm,fym,fdp,ftp,fyp=',6(1pe11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  first            = .false.
  ln_10            = DLOG(10.d+00)
END IF

!-----------------------------------------------------------------------
!        Compute independent variable grid indices and interpolation
!         coefficients.
!-----------------------------------------------------------------------

CALL eosidx_y( j, ji_ray, jk_ray, rho, t, ye, fdm, fdp, fdd, ftm, ftp, &
& ftt, fym, fyp, fyy )

!-----------------------------------------------------------------------
!        Store equation of state data for quantity i in temporary
!         scalar variables for interpolation.
!-----------------------------------------------------------------------

av111              = estble(i,j,ida  ,ita  ,iya  ,ji_ray,jk_ray)
av211              = estble(i,j,idap1,ita  ,iya  ,ji_ray,jk_ray)
av121              = estble(i,j,ida  ,itap1,iya  ,ji_ray,jk_ray)
av112              = estble(i,j,ida  ,ita  ,iyap1,ji_ray,jk_ray)
av221              = estble(i,j,idap1,itap1,iya  ,ji_ray,jk_ray)
av212              = estble(i,j,idap1,ita  ,iyap1,ji_ray,jk_ray)
av122              = estble(i,j,ida  ,itap1,iyap1,ji_ray,jk_ray)
av222              = estble(i,j,idap1,itap1,iyap1,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

besv               =  fym * ( fdm * ( ftm * av111 + ftp * av121 )   &
&                   +         fdp * ( ftm * av211 + ftp * av221 ) ) &
&                   + fyp * ( fdm * ( ftm * av112 + ftp * av122 )   &
&                   +         fdp * ( ftm * av212 + ftp * av222 ) )

IF ( DABS(besv) > alogmx ) THEN
  WRITE (nprint,1001)
  WRITE (nprint,1003)
  WRITE (nprint,1005) j, ji_ray, jk_ray, rho, t, ye
  WRITE (nprint,1007) i, besv
  WRITE (nprint,1015) av111, av211, av121, av112, av221, av212, av122, av222
  WRITE (nprint,1017) fdm, ftm, fym, fdp, ftp, fyp
  WRITE (nlog,1001)
  WRITE (nlog,1003)
  WRITE (nlog,1005) j, ji_ray, jk_ray, rho, t, ye
  WRITE (nlog,1007) i, besv
  WRITE (nlog,1015) av111, av211, av121, av112, av221, av212, av122, av222
  WRITE (nlog,1017) fdm, ftm, fym, fdp, ftp, fyp
  STOP
END IF

esvtmp             = 10.d+00**besv
esv                = esvtmp - escnst(i,j,ji_ray,jk_ray) - epsilon
IF ( i >= 7 ) esv = DMAX1( esv, zero )

!-----------------------------------------------------------------------
!  Compute d(esv)/d(density)
!-----------------------------------------------------------------------

besvd              = fdd * ( fym * ( ftm * ( -av111 + av211 )     &
&                   +                ftp * ( -av121 + av221 ) )   &
&                   +        fyp * ( ftm * ( -av112 + av212 )     &
&                   +                ftp * ( -av122 + av222 ) ) )

IF ( DABS( besvd * esvtmp ) > realmx ) THEN
  WRITE (nprint,1001)
  WRITE (nprint,1003)
  WRITE (nprint,1005) j, ji_ray, jk_ray, rho, t, ye
  WRITE (nprint,1009) i, besvd
  WRITE (nlog,1001)
  WRITE (nlog,1003)
  WRITE (nlog,1005) j, ji_ray, jk_ray, rho, t, ye
  WRITE (nlog,1009) i, besvd
  STOP
END IF

esvd               = ln_10 * besvd * esvtmp

!-----------------------------------------------------------------------
!  Compute d(esv)/d(temperature)
!-----------------------------------------------------------------------

besvt              = ftt * ( fym * ( fdm * ( -av111 + av121 )     &
&                   +                fdp * ( -av211 + av221 ) )   &
&                   +        fyp * ( fdm * ( -av112 + av122 )     &
&                   +                fdp * ( -av212 + av222 ) ) )

IF ( DABS(besvt * esvtmp) > realmx ) THEN
  WRITE (nprint,1001)
  WRITE (nprint,1003)
  WRITE (nprint,1005) j, ji_ray, jk_ray, rho, t, ye
  WRITE (nprint,1011) i, besvt
  WRITE (nlog,1001)
  WRITE (nlog,1003)
  WRITE (nlog,1005) j, ji_ray, jk_ray, rho, t, ye
  WRITE (nlog,1011) i, besvt
  STOP
END IF

esvt               = ln_10 * besvt * esvtmp

!-----------------------------------------------------------------------
!  Compute d(esv)/d(ye)
!-----------------------------------------------------------------------

besvy              = fyy * ( fdm * ( ftm * ( -av111 + av112 )     &
&                   +                ftp * ( -av121 + av122 ) )   &
&                   +        fdp * ( ftm * ( -av211 + av212 )     &
&                   +                ftp * ( -av221 + av222 ) ) )

IF ( DABS(besvy * esvtmp) > realmx ) THEN
  WRITE (nprint,1001)
  WRITE (nprint,1003)
  WRITE (nprint,1005) j, ji_ray, jk_ray, rho, t, ye
  WRITE (nprint,1013) i, besvy
  WRITE (nlog,1001)
  WRITE (nlog,1003)
  WRITE (nlog,1005) j, ji_ray, jk_ray, rho, t, ye
  WRITE (nlog,1013) i, besvy
  STOP
END IF

esvy               = ln_10 * besvy * esvtmp

RETURN
END SUBROUTINE eqstt_y

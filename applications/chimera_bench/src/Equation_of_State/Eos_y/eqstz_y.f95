SUBROUTINE eqstz_y( jmin, jmax, rho, t, ye, ji_ray, jk_ray )
!-----------------------------------------------------------------------
!
!    File:         eqstz_y
!    Module:       eqstz_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/18/05
!
!    Purpose:
!      To interpolate the logs of all equation of state quantities
!       from a local table of nearest entries created for each
!       radial shell midpoint, and to compute their derivatives
!       with respect to rho, t, and ye from the interpolation
!       formulae. This subroutine has been optimized for a vector
!       machine. This subroutine shold be called at the end of each
!       cycle of the code to update all the thermodynamic variables.
!
!    Variables that must be passed through common:
!
!  dgrid(idty(j,ji_ray,jk_ray))
!                : number of table entries per decade in rho
!                 for zone j.
!  tgrid(idty(j,ji_ray,jk_ray))
!                : number of table entries per decade in t
!                 for zone j.
!  ygrid(idty(j,ji_ray,jk_ray))
!                : number of table entries in ye between
!                 ye = 0.5 and ye = 0 for zone j.
!  idty(j,ji_ray,jk_ray) : index for dgrid, tgrid, and ygrid for zone j.
!  rhoes(i)      : density boundary between idty=i and idty=i+1
!  idr(j,ji_ray,jk_ray)  : rho grid index for zone j.
!  itr(j,ji_ray,jk_ray)  : t grid index for zone j.
!  iyr(j,ji_ray,jk_ray)  : ye grid index for zone j.
!  estble(i,j,ida,ita,iya,ji_ray,jk_ray)
!                : equation of state table array
!   i            : thermodynamic function index
!   i = 1        : pressure
!   i = 2        : energy
!   i = 3        : entropy
!   i = 4        : neutron chemical potential
!   i = 5        : proton chemical potential
!   i = 6        : electron chemical potential
!   i = 7        : free neutron mass fraction
!   i = 8        : free proton mass fraction
!   i = 9        : heavy nucleus mass fraction
!   i = 10       : heavy nucleus mass number
!   i = 11       : heavy nucleus charge number
!   i = 12       : gamma1
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jmin          : minimum y (angular) zone for which thermodynamic
!                   variables are to be evaluated
!  jmax          : maximum y (angular) zone for which thermodynamic
!                   variables are to be evaluated
!  rho           : angular array of matter density (g/cm**3)
!  t             : angular array of matter temperature (K)
!  ye            : angular array of matter electron fraction
!  ji_ray        : x (radial) index of a specific angular ray
!  jk_ray        : z (azimuthal) index of a specific angular ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : ny
USE numerical_module, ONLY : zero, one, epsilon

USE edit_module, ONLY : nprint, nlog
USE eos_snc_y_module, ONLY : idr, itr, iyr, escnst, estble, dgrid, tgrid, &
& ygrid, idty, aesv, aesvd, aesvt, aesvy

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jmin            ! minimum y (angular) zone index
INTEGER, INTENT(in)              :: jmax            ! maximum y (angular) zone index
INTEGER, INTENT(in)              :: ji_ray          ! x (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray          ! z (azimuthal) index of a specific angular ray

REAL(KIND=double), INTENT(in), DIMENSION(ny) :: rho ! angular array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: t   ! angular array of matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: ye  ! angular array of matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

LOGICAL                          :: first = .true.

INTEGER                          :: i               ! eos parameter index
INTEGER                          :: j               ! radial zone index

INTEGER                          :: id              ! density grod index
INTEGER                          :: it              ! temperature grod index
INTEGER                          :: iy              ! electron fraction grod index

INTEGER                          :: istat           ! allocation status

INTEGER, PARAMETER               :: ida = 1         ! lower cube table density index
INTEGER, PARAMETER               :: idap1 = 2       ! upper cube table density index
INTEGER, PARAMETER               :: ita = 1         ! lower cube table temperature index
INTEGER, PARAMETER               :: itap1 = 2       ! upper cube table temperature index
INTEGER, PARAMETER               :: iya = 1         ! lower cube table electron fraction index
INTEGER, PARAMETER               :: iyap1 = 2       ! upper cube table electron fraction index

REAL(KIND=double)                :: fd              ! position of rho in grid
REAL(KIND=double)                :: ft              ! position of t in grid
REAL(KIND=double)                :: fy              ! position of ye in grid

REAL(KIND=double)                :: ln_10           ! ln(10)
REAL(KIND=double)                :: loge            ! log10(e)

REAL(KIND=double)                :: av111           ! scalar table entry for interpolation
REAL(KIND=double)                :: av211           ! scalar table entry for interpolation
REAL(KIND=double)                :: av121           ! scalar table entry for interpolation
REAL(KIND=double)                :: av112           ! scalar table entry for interpolation
REAL(KIND=double)                :: av221           ! scalar table entry for interpolation
REAL(KIND=double)                :: av212           ! scalar table entry for interpolation
REAL(KIND=double)                :: av122           ! scalar table entry for interpolation
REAL(KIND=double)                :: av222           ! scalar table entry for interpolation

REAL(KIND=double)                :: besv            ! log10(esvtmp)
REAL(KIND=double)                :: besvd           ! d(besv)/d(rho)
REAL(KIND=double)                :: besvt           ! d(besv)/d(t)
REAL(KIND=double)                :: besvy           ! d(besv)/d(ye)
REAL(KIND=double)                :: esvtmp          ! esv + escnst(i,j,ji_ray,jk_ray) + epsilon

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fdp ! position of rho in grid wrt lower cube index
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fdm ! position of rho in grid wrt upper cube index
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fdd ! d(fd)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ftp ! position of t in grid wrt lower cube index
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ftm ! position of t in grid wrt upper cube index
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ftt ! d(ft)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fyp ! position of ye in grid wrt lower cube index
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fym ! position of ye in grid wrt upper cube index
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fyy ! d(fy)/d(ye)

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in eqstz_y')
 2001 FORMAT (' Deallocation problem for array ',a10,' in eqstz_y')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

  ALLOCATE (fdp(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fdm(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdm       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fdd(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdd       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftp(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftm(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftm       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftt(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftt       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fyp(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fyp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fym(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fym       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fyy(ny), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fyy       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

fdp                  = zero
fdm                  = zero
fdd                  = zero
ftp                  = zero
ftm                  = zero
ftt                  = zero
fyp                  = zero
fym                  = zero
fyy                  = zero

IF ( first ) THEN
  first              = .false.
  ln_10              = DLOG(10.d+00)
  loge               = DLOG10( DEXP(1.d+00) )
END IF

!-----------------------------------------------------------------------
!  Compute interpolation coefficients
!-----------------------------------------------------------------------

DO j = jmin,jmax

  fd                 = dgrid(idty(j,ji_ray,jk_ray)) * DLOG10(rho(j) )
  ft                 = tgrid(idty(j,ji_ray,jk_ray)) * DLOG10(t(j)   )
  fy                 = ygrid(idty(j,ji_ray,jk_ray)) * ( one - ye(j) )

  id                 = idr(j,ji_ray,jk_ray)
  fdp(j)             = fd - DBLE( id )
  fdm(j)             = one - fdp(j)
  fdd(j)             = loge * dgrid(idty(j,ji_ray,jk_ray))/( rho(j) )

  it                 = itr(j,ji_ray,jk_ray)
  ftp(j)             = ft - DBLE( it )
  ftm(j)             = one - ftp(j)
  ftt(j)             = loge * tgrid(idty(j,ji_ray,jk_ray))/( t(j)   )

  iy                 = iyr(j,ji_ray,jk_ray)
  fyp(j)             = fy - DBLE( iy )
  fym(j)             = one - fyp(j)
  fyy(j)             = - ygrid(idty(j,ji_ray,jk_ray))

END DO

!-----------------------------------------------------------------------
!
!               \\\\\ INTERPOLATE EOS PARAMETERS /////
!
!-----------------------------------------------------------------------

DO i = 1,12
  DO j = jmin,jmax

!-----------------------------------------------------------------------
!  Store equation of state data for quantity i in temporary scalar
!   variables for interpolation.
!-----------------------------------------------------------------------

    av111            = estble(i,j,ida  ,ita  ,iya  ,ji_ray,jk_ray)
    av211            = estble(i,j,idap1,ita  ,iya  ,ji_ray,jk_ray)
    av121            = estble(i,j,ida  ,itap1,iya  ,ji_ray,jk_ray)
    av112            = estble(i,j,ida  ,ita  ,iyap1,ji_ray,jk_ray)
    av221            = estble(i,j,idap1,itap1,iya  ,ji_ray,jk_ray)
    av212            = estble(i,j,idap1,ita  ,iyap1,ji_ray,jk_ray)
    av122            = estble(i,j,ida  ,itap1,iyap1,ji_ray,jk_ray)
    av222            = estble(i,j,idap1,itap1,iyap1,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Interpolate equation of state variable i (esv)
!-----------------------------------------------------------------------


    besv             =  fym(j) * ( fdm(j) * ( ftm(j) * av111 + ftp(j) * av121 )   &
&                     +            fdp(j) * ( ftm(j) * av211 + ftp(j) * av221 ) ) &
&                     + fyp(j) * ( fdm(j) * ( ftm(j) * av112 + ftp(j) * av122 )   &
&                     +            fdp(j) * ( ftm(j) * av212 + ftp(j) * av222 ) )
    esvtmp           = 10.d+00**besv
    aesv(j,i,ji_ray,jk_ray)  = esvtmp - escnst(i,j,ji_ray,jk_ray) - epsilon
    IF ( i >= 7 ) aesv(j,i,ji_ray,jk_ray) = DMAX1( aesv(j,i,ji_ray,jk_ray), zero )

!-----------------------------------------------------------------------
!  Compute d(esv)/d(density)
!-----------------------------------------------------------------------

    besvd            = fdd(j) * ( fym(j) * ( ftm(j) * ( -av111 + av211 )     &
&                     +                      ftp(j) * ( -av121 + av221 ) )   &
&                     +           fyp(j) * ( ftm(j) * ( -av112 + av212 )     &
&                     +                      ftp(j) * ( -av122 + av222 ) ) )
    aesvd(j,i,ji_ray,jk_ray) = ln_10 * besvd * esvtmp

!-----------------------------------------------------------------------
!  Compute d(esv)/d(temperature)
!-----------------------------------------------------------------------

    besvt            = ftt(j) * ( fym(j) * ( fdm(j) * ( -av111 + av121 )     &
&                     +                      fdp(j) * ( -av211 + av221 ) )   &
&                     +           fyp(j) * ( fdm(j) * ( -av112 + av122 )     &
&                     +                      fdp(j) * ( -av212 + av222 ) ) )
    aesvt(j,i,ji_ray,jk_ray) = ln_10 * besvt * esvtmp

!-----------------------------------------------------------------------
!  Compute d(esv)/d(ye)
!-----------------------------------------------------------------------

    besvy            = fyy(j) * ( fdm(j) * ( ftm(j) * ( -av111 + av112 )     &
&                     +                      ftp(j) * ( -av121 + av122 ) )   &
&                     +           fdp(j) * ( ftm(j) * ( -av211 + av212 )     &
&                     +                      ftp(j) * ( -av221 + av222 ) ) )
    aesvy(j,i,ji_ray,jk_ray) = ln_10 * besvy * esvtmp

  END DO
END DO

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

  DEALLOCATE (fdp, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdp       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (fdm, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdm       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (fdd, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdd       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (ftp, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftp       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (ftm, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftm       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (ftt, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftt       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (fyp, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fyp       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (fym, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fym       '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (fyy, STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fyy       '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE eqstz_y

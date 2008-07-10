SUBROUTINE eqstz_z( kmin, kmax, rho, t, ye, ki_ray, kj_ray )
!-----------------------------------------------------------------------
!
!    File:         eqstz_z
!    Module:       eqstz_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
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
!  dgrid(idty(k,kj_ray,ki_ray))
!                : number of table entries per decade in rho
!                 for zone k.
!  tgrid(idty(k,kj_ray,ki_ray))
!                : number of table entries per decade in t
!                 for zone k.
!  ygrid(idty(k,kj_ray,ki_ray))
!                : number of table entries in ye between
!                 ye = 0.5 and ye = 0 for zone k.
!  idty(k,kj_ray,ki_ray) : index for dgrid, tgrid, and ygrid for zone k.
!  rhoes(i)      : density boundary between idty=i and idty=i+1
!  idr(k,kj_ray,ki_ray)  : rho grid index for zone k.
!  itr(k,kj_ray,ki_ray)  : t grid index for zone k.
!  iyr(k,kj_ray,ki_ray)  : ye grid index for zone k.
!  estble(i,k,ida,ita,iya,kj_ray,ki_ray)
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
!  kmin          : minimum z (azimuthal) zone for which thermodynamic
!                   variables are to be evaluated
!  kmax          : maximum z (azimuthal) zone for which thermodynamic
!                   variables are to be evaluated
!  rho           : azimuthal array of matter density (g/cm**3)
!  t             : azimuthal array of matter temperature (K)
!  ye            : azimuthal array of matter electron fraction
!  ki_ray        : x (radial) index of a specific z (azimuthaal) ray
!  kj_ray        : y (angular) index of a specific z (azimuthaal) ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_z_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nz
USE numerical_module, ONLY : zero, one, epsilon

USE edit_module, ONLY : nprint, nlog
USE eos_snc_z_module, ONLY : idr, itr, iyr, escnst, estble, dgrid, tgrid, &
& ygrid, idty, aesv, aesvd, aesvt, aesvy

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: kmin            ! minimum z (azimuthal) zone index
INTEGER, INTENT(in)              :: kmax            ! maximum z (azimuthal) zone index
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthaal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthaal) ray

REAL(KIND=double), INTENT(in), DIMENSION(nz) :: rho ! azimuthal array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: t   ! azimuthal array of matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: ye  ! azimuthal array of matter matter electron fraction array

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

LOGICAL                          :: first = .true.

INTEGER                          :: i               ! eos parameter index
INTEGER                          :: k               ! z (azimuthal) zone index.

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
REAL(KIND=double)                :: esvtmp          ! esv + escnst(i,k,kj_ray,ki_ray) + epsilon

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

 1001 FORMAT (' Allocation problem for array ',a10,' in eqstz_z')
 2001 FORMAT (' Deallocation problem for array ',a10,' in eqstz_z')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

  ALLOCATE (fdp(nz), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fdm(nz), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdm       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fdd(nz), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdd       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftp(nz), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftm(nz), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftm       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftt(nz), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftt       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fyp(nz), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fyp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fym(nz), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fym       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fyy(nz), STAT = istat)
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

DO k = kmin,kmax

  fd                 = dgrid(idty(k,kj_ray,ki_ray)) * DLOG10(rho(k) )
  ft                 = tgrid(idty(k,kj_ray,ki_ray)) * DLOG10(t(k)   )
  fy                 = ygrid(idty(k,kj_ray,ki_ray)) * ( one - ye(k) )

  id                 = idr(k,kj_ray,ki_ray)
  fdp(k)             = fd - DBLE( id )
  fdm(k)             = one - fdp(k)
  fdd(k)             = loge * dgrid(idty(k,kj_ray,ki_ray))/( rho(k) )

  it                 = itr(k,kj_ray,ki_ray)
  ftp(k)             = ft - DBLE( it )
  ftm(k)             = one - ftp(k)
  ftt(k)             = loge * tgrid(idty(k,kj_ray,ki_ray))/( t(k)   )

  iy                 = iyr(k,kj_ray,ki_ray)
  fyp(k)             = fy - DBLE( iy )
  fym(k)             = one - fyp(k)
  fyy(k)             = - ygrid(idty(k,kj_ray,ki_ray))

END DO

!-----------------------------------------------------------------------
!
!               \\\\\ INTERPOLATE EOS PARAMETERS /////
!
!-----------------------------------------------------------------------

DO i = 1,12
  DO k = kmin,kmax

!-----------------------------------------------------------------------
!  Store equation of state data for quantity i in temporary scalar
!   variables for interpolation.
!-----------------------------------------------------------------------

    av111            = estble(i,k,ida  ,ita  ,iya  ,kj_ray,ki_ray)
    av211            = estble(i,k,idap1,ita  ,iya  ,kj_ray,ki_ray)
    av121            = estble(i,k,ida  ,itap1,iya  ,kj_ray,ki_ray)
    av112            = estble(i,k,ida  ,ita  ,iyap1,kj_ray,ki_ray)
    av221            = estble(i,k,idap1,itap1,iya  ,kj_ray,ki_ray)
    av212            = estble(i,k,idap1,ita  ,iyap1,kj_ray,ki_ray)
    av122            = estble(i,k,ida  ,itap1,iyap1,kj_ray,ki_ray)
    av222            = estble(i,k,idap1,itap1,iyap1,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Interpolate equation of state variable i (esv)
!-----------------------------------------------------------------------


    besv             =  fym(k) * ( fdm(k) * ( ftm(k) * av111 + ftp(k) * av121 )   &
&                     +            fdp(k) * ( ftm(k) * av211 + ftp(k) * av221 ) ) &
&                     + fyp(k) * ( fdm(k) * ( ftm(k) * av112 + ftp(k) * av122 )   &
&                     +            fdp(k) * ( ftm(k) * av212 + ftp(k) * av222 ) )
    esvtmp           = 10.d+00**besv
    aesv(k,i,kj_ray,ki_ray)  = esvtmp - escnst(i,k,kj_ray,ki_ray) - epsilon
    IF ( i >= 7 ) aesv(k,i,kj_ray,ki_ray) = DMAX1( aesv(k,i,kj_ray,ki_ray), zero )

!-----------------------------------------------------------------------
!  Compute d(esv)/d(density)
!-----------------------------------------------------------------------

    besvd            = fdd(k) * ( fym(k) * ( ftm(k) * ( -av111 + av211 )     &
&                     +                      ftp(k) * ( -av121 + av221 ) )   &
&                     +           fyp(k) * ( ftm(k) * ( -av112 + av212 )     &
&                     +                      ftp(k) * ( -av122 + av222 ) ) )
    aesvd(k,i,kj_ray,ki_ray) = ln_10 * besvd * esvtmp

!-----------------------------------------------------------------------
!  Compute d(esv)/d(temperature)
!-----------------------------------------------------------------------

    besvt            = ftt(k) * ( fym(k) * ( fdm(k) * ( -av111 + av121 )     &
&                     +                      fdp(k) * ( -av211 + av221 ) )   &
&                     +           fyp(k) * ( fdm(k) * ( -av112 + av122 )     &
&                     +                      fdp(k) * ( -av212 + av222 ) ) )
    aesvt(k,i,kj_ray,ki_ray) = ln_10 * besvt * esvtmp

!-----------------------------------------------------------------------
!  Compute d(esv)/d(ye)
!-----------------------------------------------------------------------

    besvy            = fyy(k) * ( fdm(k) * ( ftm(k) * ( -av111 + av112 )     &
&                     +                      ftp(k) * ( -av121 + av122 ) )   &
&                     +           fdp(k) * ( ftm(k) * ( -av211 + av212 )     &
&                     +                      ftp(k) * ( -av221 + av222 ) ) )
    aesvy(k,i,kj_ray,ki_ray) = ln_10 * besvy * esvtmp

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
END SUBROUTINE eqstz_z

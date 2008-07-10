SUBROUTINE eqstz_x( jr_min, jr_max, rho, t, ye, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         eqstz_x
!    Module:       eqstz_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
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
!  dgrid(idty(j,i_ray))
!                : number of table entries per decade in rho
!                   for zone j.
!  tgrid(idty(j,i_ray))
!                : number of table entries per decade in t
!                   for zone j.
!  ygrid(idty(j,i_ray))
!                : number of table entries in ye between
!                   ye = 0.5 and ye = 0 for zone j.
!  idty(j,i_ray) : index for dgrid, tgrid, and ygrid for zone j.
!  rhoes(i)      : density boundary between idty=i and idty=i+1
!  idr(j,i_ray)  : rho grid index for zone j.
!  itr(j,i_ray)  : t grid index for zone j.
!  iyr(j,i_ray)  : ye grid index for zone j.
!  estble(i,j,ida,ita,iya,i_ray)
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
!  jr_min        : minimum radial zone for which thermodynamic
!                   variables are to be evaluated
!  jr_max        : maximum radial zone for which thermodynamic
!                   variables are to be evaluated
!  rho           : shifted matter density array (g/cm**3).
!  t             : shifted matter matter temperature array (K).
!  ye            : shifted matter matter electron fraction array.
!  ij_ray        : j-index of a radial ray
!  ik_ray        : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx
USE numerical_module, ONLY : zero, one, epsilon

USE edit_module, ONLY : nprint, nlog
USE eos_snc_x_module, ONLY : idr, itr, iyr, escnst, estble, dgrid, tgrid, &
& ygrid, idty, aesv, aesvd, aesvt, aesvy

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min          ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max          ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t   ! shifted matter matter temperature array (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: ye  ! shifted matter matter electron fraction array

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
REAL(KIND=double)                :: esvtmp          ! esv + escnst(i,j,ij_ray,ik_ray) + epsilon

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

 1001 FORMAT (' Allocation problem for array ',a10,' in eqstz_x')
 2001 FORMAT (' Deallocation problem for array ',a10,' in eqstz_x')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

  ALLOCATE (fdp(nx), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fdm(nx), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdm       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fdd(nx), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fdd       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftp(nx), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftm(nx), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftm       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (ftt(nx), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'ftt       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fyp(nx), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fyp       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fym(nx), STAT = istat)
   IF ( istat /= 0 ) THEN; var_name = 'fym       '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (fyy(nx), STAT = istat)
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

DO j = jr_min,jr_max

  fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10(rho(j) )
  ft                 = tgrid(idty(j,ij_ray,ik_ray)) * DLOG10(t(j)   )
  fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye(j) )

  id                 = idr(j,ij_ray,ik_ray)
  fdp(j)             = fd - DBLE( id )
  fdm(j)             = one - fdp(j)
  fdd(j)             = loge * dgrid(idty(j,ij_ray,ik_ray))/( rho(j) )

  it                 = itr(j,ij_ray,ik_ray)
  ftp(j)             = ft - DBLE( it )
  ftm(j)             = one - ftp(j)
  ftt(j)             = loge * tgrid(idty(j,ij_ray,ik_ray))/( t(j)   )

  iy                 = iyr(j,ij_ray,ik_ray)
  fyp(j)             = fy - DBLE( iy )
  fym(j)             = one - fyp(j)
  fyy(j)             = - ygrid(idty(j,ij_ray,ik_ray))

END DO

!-----------------------------------------------------------------------
!
!               \\\\\ INTERPOLATE EOS PARAMETERS /////
!
!-----------------------------------------------------------------------

DO i = 1,12
  DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Store equation of state data for quantity i in temporary scalar
!   variables for interpolation.
!-----------------------------------------------------------------------

    av111            = estble(i,j,ida  ,ita  ,iya  ,ij_ray,ik_ray)
    av211            = estble(i,j,idap1,ita  ,iya  ,ij_ray,ik_ray)
    av121            = estble(i,j,ida  ,itap1,iya  ,ij_ray,ik_ray)
    av112            = estble(i,j,ida  ,ita  ,iyap1,ij_ray,ik_ray)
    av221            = estble(i,j,idap1,itap1,iya  ,ij_ray,ik_ray)
    av212            = estble(i,j,idap1,ita  ,iyap1,ij_ray,ik_ray)
    av122            = estble(i,j,ida  ,itap1,iyap1,ij_ray,ik_ray)
    av222            = estble(i,j,idap1,itap1,iyap1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Interpolate equation of state variable i (esv)
!-----------------------------------------------------------------------


    besv             =  fym(j) * ( fdm(j) * ( ftm(j) * av111 + ftp(j) * av121 )   &
&                     +            fdp(j) * ( ftm(j) * av211 + ftp(j) * av221 ) ) &
&                     + fyp(j) * ( fdm(j) * ( ftm(j) * av112 + ftp(j) * av122 )   &
&                     +            fdp(j) * ( ftm(j) * av212 + ftp(j) * av222 ) )
    esvtmp           = 10.d+00**besv
    aesv(j,i,ij_ray,ik_ray)  = esvtmp - escnst(i,j,ij_ray,ik_ray) - epsilon
    IF ( i >= 7 ) aesv(j,i,ij_ray,ik_ray) = DMAX1( aesv(j,i,ij_ray,ik_ray), zero )

!-----------------------------------------------------------------------
!  Compute d(esv)/d(density)
!-----------------------------------------------------------------------

    besvd            = fdd(j) * ( fym(j) * ( ftm(j) * ( -av111 + av211 )     &
&                     +                      ftp(j) * ( -av121 + av221 ) )   &
&                     +           fyp(j) * ( ftm(j) * ( -av112 + av212 )     &
&                     +                      ftp(j) * ( -av122 + av222 ) ) )
    aesvd(j,i,ij_ray,ik_ray) = ln_10 * besvd * esvtmp

!-----------------------------------------------------------------------
!  Compute d(esv)/d(temperature)
!-----------------------------------------------------------------------

    besvt            = ftt(j) * ( fym(j) * ( fdm(j) * ( -av111 + av121 )     &
&                     +                      fdp(j) * ( -av211 + av221 ) )   &
&                     +           fyp(j) * ( fdm(j) * ( -av112 + av122 )     &
&                     +                      fdp(j) * ( -av212 + av222 ) ) )
    aesvt(j,i,ij_ray,ik_ray) = ln_10 * besvt * esvtmp

!-----------------------------------------------------------------------
!  Compute d(esv)/d(ye)
!-----------------------------------------------------------------------

    besvy            = fyy(j) * ( fdm(j) * ( ftm(j) * ( -av111 + av112 )     &
&                     +                      ftp(j) * ( -av121 + av122 ) )   &
&                     +           fdp(j) * ( ftm(j) * ( -av211 + av212 )     &
&                     +                      ftp(j) * ( -av221 + av222 ) ) )
    aesvy(j,i,ij_ray,ik_ray) = ln_10 * besvy * esvtmp

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
END SUBROUTINE eqstz_x

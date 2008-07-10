SUBROUTINE bremkrnl( n, j, ij_ray, ik_ray, k , rho, t, ye, f0a, f1a, &
& f0ad, f1ad, f0at, f1at, f0ay, f1ay, nez )
!-----------------------------------------------------------------------
!
!    File:         bremkrnl
!    Module:       bremkrnl
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/6/02
!
!    Purpose:
!      To interpolate pair annihilation bremsstrahlung kernals from a local
!       table of nearest entries  created for each zone.
!      The table consists of the eight nearest-neighbor entries in rho, t,
!       and ye of the state point.
!      Derivatives of the pair annihilation rates are obtained from the
!       interpolation formula for these rates by direct differentiation
!       of the interpolation formula.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  n                    : neutrino type
!  j                    : radial zone number
!  ij_ray               : j-index of a radial ray
!  ik_ray               : k-index of a radial ray
!  k                    : neutrino energy group
!  rho                  : density (g/cm**3)
!  t                    : temperature (K)
!  ye                   : electron fraction
!  nez                  : neutrino energy array dimension
!
!    Output arguments:
!  f0a                  : zero moment of the bremsstrahlung pair annihilation kernal
!  f01                  : first moment of the bremsstrahlung pair annihilation kernal
!  f0ad                 : derivative with respect to the density of the zero
!                          moment of the bremsstrahlung pair annihilation kernal
!  f1ad                 : derivative with respect to the density of the first moment
!                          of the bremsstrahlung pair annihilation kernal
!  f0at                 : derivative with respect to the temperature of the zero moment
!                          of the bremsstrahlung pair annihilation kernal
!  f1at                 : derivative with respect to the temperature of the first moment
!                          of the bremsstrahlung pair annihilation kernal
!  f0ay                 : derivative with respect to the electron fraction of the zero
!                          moment of the bremsstrahlung pair annihilation kernal
!  f1ay                 : derivative with respect to the electron fraction of the first
!                          moment of the pair bremsstrahlung annihilation kernal
!
!    Variables that must be passed through common:
!
!  rhobrememn           : rho < rhobrememn: e-neutrino-antineutrino bremsstrahlung pair
!                          annihilation and production omitted
!  rhobrememx           : rho > rhobrememx: e-neutrino-antineutrino bremsstrahlung pair
!                          annihilation and production omitted
!  rhobremtmn           : rho < rhobremtmn: t-neutrino-antineutrino pair bremsstrahlung
!                          pair annihilation and production omitted
!  rhobremtmx           : rho > rhobremtmx: t-neutrino-antineutrino pair bremsstrahlung
!                          pair annihilation and production omitted
!  brema0               : zero moment of the neutrino-antineutrino bremsstrahlung
!                          annihilation kernal array
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  idrb(j,ij_ray,ik_ray)        : rho grid index for zone j
!  itrb(j,ij_ray,ik_ray)        : t grid index for zone j
!  iyrb(j,ij_ray,ik_ray)        : ye grid index for zone j
!
!    Include files:
!  kind_module, array_module, numerical_module
!  brem_module, eos_snc_x_module, nu_energy_grid_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : ij_ray_dim
USE numerical_module, ONLY : zero, half, one

USE brem_module, ONLY : brema0, idrb, itrb, iyrb
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY : rhobrememn, rhobrememx, rhobremtmn, rhobremtmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: j             ! radial zone index
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)               :: k             ! incoming neutrino energy zone index
INTEGER, INTENT(in)               :: nez           ! neutrino energy array dimension

REAL(KIND=double), INTENT(in)     :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)     :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)     :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0a  ! zero moment of the pair annihilation function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f1a  ! first moment of the pair annihilation function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0ad ! d(f0a)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f1ad ! d(f1a)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0at ! d(f0a)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f1at ! d(f1a)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0ay ! d(f0a)/d(ye)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f1ay ! d(f1a)/d(ye)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.
LOGICAL                          :: brem_off

INTEGER                          :: kp            ! antineutrino energy zone index

INTEGER                          :: id            ! EOS density index
INTEGER                          :: it            ! EOS temperature index
INTEGER                          :: iy            ! EOS electron fraction index

INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

INTEGER, PARAMETER               :: ida = 1       ! lower cube table density index
INTEGER, PARAMETER               :: idap1 = 2     ! upper cube table density index
INTEGER, PARAMETER               :: ita = 1       ! lower cube table temperature index
INTEGER, PARAMETER               :: itap1 = 2     ! upper cube table temperature index
INTEGER, PARAMETER               :: iya = 1       ! lower cube table electron fraction index
INTEGER, PARAMETER               :: iyap1 = 2     ! upper cube table electron fraction index

REAL(KIND=double)                :: log_e         ! log10(e)
REAL(KIND=double)                :: ln_10         ! ln(10)

REAL(KIND=double)                :: fd            ! position of rho in grid
REAL(KIND=double)                :: fdp           ! position of rho in grid wrt lower cube index
REAL(KIND=double)                :: fdm           ! position of rho in grid wrt upper cube index
REAL(KIND=double)                :: fdd           ! d(fd)/d(rho)
REAL(KIND=double)                :: ft            ! position of t in grid
REAL(KIND=double)                :: ftp           ! position of t in grid wrt lower cube index
REAL(KIND=double)                :: ftm           ! position of t in grid wrt upper cube index
REAL(KIND=double)                :: ftt           ! d(ft)/d(t)
REAL(KIND=double)                :: fy            ! position of ye in grid
REAL(KIND=double)                :: fyp           ! position of ye in grid wrt lower cube index
REAL(KIND=double)                :: fym           ! position of ye in grid wrt upper cube index
REAL(KIND=double)                :: fyy           ! d(fy)/d(ye)

REAL(KIND=double)                :: b111          ! scalar table entry for interpolation
REAL(KIND=double)                :: b211          ! scalar table entry for interpolation
REAL(KIND=double)                :: b121          ! scalar table entry for interpolation
REAL(KIND=double)                :: b112          ! scalar table entry for interpolation
REAL(KIND=double)                :: b221          ! scalar table entry for interpolation
REAL(KIND=double)                :: b212          ! scalar table entry for interpolation
REAL(KIND=double)                :: b122          ! scalar table entry for interpolation
REAL(KIND=double)                :: b222          ! scalar table entry for interpolation

REAL(KIND=double)                :: f0al          ! interpolated log of zero-m B annihilation funciton
REAL(KIND=double)                :: f0adl         ! interpolated log of d(f0a)/d(rho)
REAL(KIND=double)                :: f0atl         ! interpolated log of d(f0a)/d(t)
REAL(KIND=double)                :: f0ayl         ! interpolated log of d(f0a)/d(ye)

REAL(KIND=double)                :: f10           ! 10** function
EXTERNAL f10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) RETURN

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF (first) THEN
  first            = .false.
  log_e            = DLOG10( DEXP( one ) )
  ln_10            = DLOG( 1.d+01 )
END IF

!-----------------------------------------------------------------------
!  Set bremsstrahlung pair annihilation functions to zero if specified
!   boundaries.
!-----------------------------------------------------------------------

brem_off           = .false.

IF ( n < 3   .and.  rho < rhobrememn ) brem_off = .true.
IF ( n < 3   .and.  rho > rhobrememx ) brem_off = .true.
IF ( n >= 3  .and.  rho < rhobremtmn ) brem_off = .true.
IF ( n >= 3  .and.  rho > rhobremtmx ) brem_off = .true.
IF ( brem_off ) THEN

  DO kp = 1,nnugp(n)
    f0a (kp)       = zero
    f0ad(kp)       = zero
    f0at(kp)       = zero
    f0ay(kp)       = zero
    f1a (kp)       = zero
    f1ad(kp)       = zero
    f1at(kp)       = zero
    f1ay(kp)       = zero
  END DO
RETURN

END IF ! brem_off

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho )
id                 = idrb(j,ij_ray,ik_ray)
fdp                = DMAX1( DMIN1( fd - DBLE( id ), one ), zero )
fdm                = one - fdp
fdd                = log_e * dgrid(idty(j,ij_ray,ik_ray))/( rho ) ! log_e = log(e)

ft                 = tgrid(idty(j,ij_ray,ik_ray)) * DLOG10( t   )
it                 = itrb(j,ij_ray,ik_ray)
ftp                = DMAX1( DMIN1( ft - DBLE( it ), one ), zero )
ftm                = one - ftp
ftt                = log_e * tgrid(idty(j,ij_ray,ik_ray))/( t   ) ! log_e = log(e)

fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye )
iy                 = iyrb(j,ij_ray,ik_ray)
fyp                = DMAX1( DMIN1( fy - DBLE( iy ), one ), zero )
fym                = one - fyp
fyy                = - ygrid(idty(j,ij_ray,ik_ray))

DO kp = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store neutrino-antineutrino pair annihilation functions in temporary
!   scalar variables for interpolation
!-----------------------------------------------------------------------

  b111             = brema0(j,k,kp,ida  ,ita  ,iya  ,i_ray)
  b211             = brema0(j,k,kp,idap1,ita  ,iya  ,i_ray)
  b121             = brema0(j,k,kp,ida  ,itap1,iya  ,i_ray)
  b112             = brema0(j,k,kp,ida  ,ita  ,iyap1,i_ray)
  b221             = brema0(j,k,kp,idap1,itap1,iya  ,i_ray)
  b212             = brema0(j,k,kp,idap1,ita  ,iyap1,i_ray)
  b122             = brema0(j,k,kp,ida  ,itap1,iyap1,i_ray)
  b222             = brema0(j,k,kp,idap1,itap1,iyap1,i_ray)

!-----------------------------------------------------------------------
!  Interpolate neutrino-antineutrino bremsstrahlung pair annihilation
!   kernals.
!-----------------------------------------------------------------------

  f0al             = fym * ( fdm * ( ftm * b111  + ftp * b121  )   &
&                  +         fdp * ( ftm * b211  + ftp * b221  ) ) &
&                  + fyp * ( fdm * ( ftm * b112  + ftp * b122  )   &
&                  +         fdp * ( ftm * b212  + ftp * b222  ) )

  f0a (kp)           = f10(f0al)
  f1a (kp)           = zero

!-----------------------------------------------------------------------
!  Interpolate d(bremsstrahlung pair annihilation kernals)/d(density)
!-----------------------------------------------------------------------

  f0adl            = fdd * ( fym * ( ftm * ( -b111  + b211  )    &
&                  +                 ftp * ( -b121  + b221  ) )  &
&                  +         fyp * ( ftm * ( -b112  + b212  )    &
&                  +                 ftp * ( -b122  + b222  ) ) )

  f0ad(kp)         = ln_10 * f0adl * f0a(kp)
  f1ad(kp)         = zero

!-----------------------------------------------------------------------
!  Interpolate d(bremsstrahlung pair annihilation kernals)/d(temperature)
!-----------------------------------------------------------------------
  f0atl            = ftt * ( fym * ( fdm * ( -b111  + b121  )    &
&                  +                 fdp * ( -b211  + b221  ) )  &
&                  +         fyp * ( fdm * ( -b112  + b122  )    &
&                  +                 fdp * ( -b212  + b222  ) ) )

  f0at(kp)         = ln_10 * f0atl * f0a(kp)
  f1at(kp)         = zero

!-----------------------------------------------------------------------
!  Interpolate d(bremsstrahlung pair annihilation kernals)/d(ye)
!-----------------------------------------------------------------------

  f0ayl            = fyy * ( fdm * ( ftm * ( -b111  + b112  )    &
&                  +                 ftp * ( -b121  + b122  ) )  &
&                  +         fdp * ( ftm * ( -b211  + b212  )    &
&                  +                 ftp * ( -b221  + b222  ) ) )

  f0ay(kp)         = ln_10 * f0ayl * f0a(kp)
  f1ay(kp)         = zero

END DO

RETURN
END SUBROUTINE bremkrnl

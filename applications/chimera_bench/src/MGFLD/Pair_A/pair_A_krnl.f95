SUBROUTINE pair_A_krnl( n, j, ij_ray, ik_ray, k, rho, t, ye, f0a, f0p, &
& f0ad, f0pd, f0at, f0pt, f0ay, f0py, nez )
!-----------------------------------------------------------------------
!
!    File:         pair_A_krnl
!    Module:       pair_A_krnl
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/4/03
!
!    Purpose:
!      To interpolate pair annihilation scattering kernals from a
!       local table of nearest entries created for each zone.
!      The table consists of the eight nearest-neighbor entries in
!       rho, t, and ye of the state point.
!      Derivatives of the pair annihilation rates are obtained
!       from the interpolation formula for these rates by direct
!       differentiation of the interpolation formula.
!
!    Input arguments:
!  n               : neutrino type
!  j               : radial zone number
!  ij_ray          : j-index of a radial ray
!  ik_ray          : k-index of a radial ray
!  k               : neutrino energy grouip
!  nez             : neutrino energy array dimension
!  rho             : density (g/cm**3)
!  t               : temperature (K)
!  ye              : electron fraction
!
!    Output arguments:
!  f0a             : zero moment of the pair annihilation kernal
!  f01             : first moment of the pair annihilation kernal
!  f0ad            : derivative with respect to the density of the
!                     zero moment of the pair annihilation kernal
!  f0pd            : derivative with respect to the density of the
!                     first moment of the pair annihilation kernal
!  f0at            : derivative with respect to the temperature of the
!                     zero moment of the pair annihilation kernal
!  f0pt            : derivative with respect to the temperature of the
!                     first moment of the pair annihilation kernal
!  f0ay            : derivative with respect to the electron fraction
!                     of the zero moment of the pair annihilation kernal
!  f0py            : derivative with respect to the electron fraction
!                     of the first moment of the pair annihilation kernal
!
!    Variables that must be passed through common:
!
!  rhopairAemn     : rho < rhopairAemn: e-neutrino-antineutrino pair annihilation and production omitted
!  rhopairAemx     : rho > rhopairAemx: e-neutrino-antineutrino pair annihilation and production omitted
!  rhopairAxmn     : rho < rhopairAxmn: x-neutrino-antineutrino pair annihilation and production omitted
!  rhopairAxmx     : rho > rhopairAxmx: x-neutrino-antineutrino pair annihilation and production omitted
!  pair_A_a        : zero moment of the nuclar neutrino-antineutrino pair annihilation kernal array
!  pair_A_p        : zero moment of the nuclar neutrino-antineutrino pair production kernal array
!  dgrid(idty(j))  : number of table entries per decade in rho for zone j
!  tgrid(idty(j))  : number of table entries per decade in t for zone j
!  ygrid(idty(j))  : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)   : index for dgrid, tgrid, and ygrid for zone j
!  idrp_A(j,ij_ray,ik_ray) : rho grid index for zone j
!  itrp_A(j,ij_ray,ik_ray) : t grid index for zone j
!  iyrp_A(j,ij_ray,ik_ray) : ye grid index for zone j
!
!    Subprograms called:
!      none
!
!    Include files:
!      array_module, kind_module, numerical_module
!      eos_snc_x_module, nu_energy_grid_module, pair_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu, ij_ray_dim
USE numerical_module, ONLY : zero, half, one

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp
USE pair_A_module, ONLY : idrp_A, itrp_A, iyrp_A, pair_A_a, pair_A_p
USE prb_cntl_module, ONLY : rhopairAemn, rhopairAemx, rhopairAtmn, rhopairAtmx

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

REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0a  ! zero moment of the nuclear pair annihilation function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0p  ! zero moment of the nuclear pair production function
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0ad ! d(f0a)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0pd ! d(f0p)/d(rho)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0at ! d(f0a)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0pt ! d(f0p)/d(t)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0ay ! d(f0a)/d(ye)
REAL(KIND=double), INTENT(out), DIMENSION(nez) :: f0py ! d(f0p)/d(ye)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.
LOGICAL                          :: pair_off

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

REAL(KIND=double), DIMENSION(nnu) :: cpair1       ! neutrino flavor coefficient
REAL(KIND=double), DIMENSION(nnu) :: cpair2       ! neutrino flavor coefficient

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

REAL(KIND=double)                :: pa111         ! scalar table entry for interpolation
REAL(KIND=double)                :: pa211         ! scalar table entry for interpolation
REAL(KIND=double)                :: pa121         ! scalar table entry for interpolation
REAL(KIND=double)                :: pa112         ! scalar table entry for interpolation
REAL(KIND=double)                :: pa221         ! scalar table entry for interpolation
REAL(KIND=double)                :: pa212         ! scalar table entry for interpolation
REAL(KIND=double)                :: pa122         ! scalar table entry for interpolation
REAL(KIND=double)                :: pa222         ! scalar table entry for interpolation

REAL(KIND=double)                :: pp111         ! scalar table entry for interpolation
REAL(KIND=double)                :: pp211         ! scalar table entry for interpolation
REAL(KIND=double)                :: pp121         ! scalar table entry for interpolation
REAL(KIND=double)                :: pp112         ! scalar table entry for interpolation
REAL(KIND=double)                :: pp221         ! scalar table entry for interpolation
REAL(KIND=double)                :: pp212         ! scalar table entry for interpolation
REAL(KIND=double)                :: pp122         ! scalar table entry for interpolation
REAL(KIND=double)                :: pp222         ! scalar table entry for interpolation

REAL(KIND=double)                :: fa0l          ! interpolated log of zero-m pair funciton i
REAL(KIND=double)                :: fp0l          ! interpolated log of zero-m pair funciton i
REAL(KIND=double)                :: fa0           ! interpolated zero-m pair funciton i
REAL(KIND=double)                :: fp0           ! interpolated zero-m pair funciton i

REAL(KIND=double)                :: fa0dl         ! interpolated log of d(fa0)/d(rho)
REAL(KIND=double)                :: fp0dl        ! interpolated log of d(fp0)/d(rho)
REAL(KIND=double)                :: fa0d          ! interpolated d(fa0)/d(rho)
REAL(KIND=double)                :: fp0d         ! interpolated d(fp0)/d(rho)

REAL(KIND=double)                :: fa0tl         ! interpolated log of d(fa0)/d(t)
REAL(KIND=double)                :: fp0tl        ! interpolated log of d(fp0)/d(t)
REAL(KIND=double)                :: fa0t          ! interpolated d(fa0)/d(t)
REAL(KIND=double)                :: fp0t         ! interpolated d(jii)/d(t)

REAL(KIND=double)                :: fa0yl         ! interpolated log of d(fa0)/d(ye)
REAL(KIND=double)                :: fp0yl        ! interpolated log of d(jii)/d(ye)
REAL(KIND=double)                :: fa0y          ! interpolated d(fa0)/d(ye)
REAL(KIND=double)                :: fp0y         ! interpolated d(fp0)/d(ye)

REAL(KIND=double)                :: f10           ! 10** function

EXTERNAL f10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) .eq. 0 ) RETURN

!-----------------------------------------------------------------------
!  Compute standard model coefficients of the neutrino-nucleus pair
!   annihilation functions.
!-----------------------------------------------------------------------

IF (first) THEN
  first            = .false.
  log_e            = DLOG10( DEXP( one ) )
  ln_10            = DLOG( 1.d+01 )
END IF

!-----------------------------------------------------------------------
!  Set pair annihilation functions to zero if rho outside specified
!   boundaries.
!-----------------------------------------------------------------------

pair_off           = .false.

IF ( n < 3  .and.  rho < rhopairAemn ) pair_off = .true.
IF ( n < 3  .and.  rho > rhopairAemx ) pair_off = .true.
IF ( n == 3  .and.  rho < rhopairAtmn ) pair_off = .true.
IF ( n == 3  .and.  rho > rhopairAtmx ) pair_off = .true.

IF ( pair_off ) THEN

  DO kp = 1,nnugp(n)
    f0a (kp)       = zero
    f0ad(kp)       = zero
    f0at(kp)       = zero
    f0ay(kp)       = zero
    f0p (kp)       = zero
    f0pd(kp)       = zero
    f0pt(kp)       = zero
    f0py(kp)       = zero
  END DO
  RETURN

END IF ! pair_off

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho )
id                 = idrp_A(j,ij_ray,ik_ray)
fdp                = DMAX1( DMIN1( fd - DBLE( id ), one ), zero )
fdm                = one - fdp
fdd                = log_e * dgrid(idty(j,ij_ray,ik_ray))/( rho ) ! log_e = log(e)

ft                 = tgrid(idty(j,ij_ray,ik_ray)) * dlog10( t   )
it                 = itrp_A(j,ij_ray,ik_ray)
ftp                = DMAX1( DMIN1( ft - DBLE( it ), one ), zero )
ftm                = one - ftp
ftt                = log_e * tgrid(idty(j,ij_ray,ik_ray))/( t   ) ! log_e = log(e)

fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye )
iy                 = iyrp_A(j,ij_ray,ik_ray)
fyp                = DMAX1( DMIN1( fy - DBLE( iy ), one ), zero )
fym                = one - fyp
fyy                = - ygrid(idty(j,ij_ray,ik_ray))

DO kp = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store neutrino-nucleus pair annihilation functions in temporary
!   scalar variables for interpolation
!-----------------------------------------------------------------------

  pa111            = pair_A_a(j,k,kp,ida  ,ita  ,iya  ,i_ray)
  pa211            = pair_A_a(j,k,kp,idap1,ita  ,iya  ,i_ray)
  pa121            = pair_A_a(j,k,kp,ida  ,itap1,iya  ,i_ray)
  pa112            = pair_A_a(j,k,kp,ida  ,ita  ,iyap1,i_ray)
  pa221            = pair_A_a(j,k,kp,idap1,itap1,iya  ,i_ray)
  pa212            = pair_A_a(j,k,kp,idap1,ita  ,iyap1,i_ray)
  pa122            = pair_A_a(j,k,kp,ida  ,itap1,iyap1,i_ray)
  pa222            = pair_A_a(j,k,kp,idap1,itap1,iyap1,i_ray)

  pp111            = pair_A_p(j,k,kp,ida  ,ita  ,iya  ,i_ray)
  pp211            = pair_A_p(j,k,kp,idap1,ita  ,iya  ,i_ray)
  pp121            = pair_A_p(j,k,kp,ida  ,itap1,iya  ,i_ray)
  pp112            = pair_A_p(j,k,kp,ida  ,ita  ,iyap1,i_ray)
  pp221            = pair_A_p(j,k,kp,idap1,itap1,iya  ,i_ray)
  pp212            = pair_A_p(j,k,kp,idap1,ita  ,iyap1,i_ray)
  pp122            = pair_A_p(j,k,kp,ida  ,itap1,iyap1,i_ray)
  pp222            = pair_A_p(j,k,kp,idap1,itap1,iyap1,i_ray)

!-----------------------------------------------------------------------
!  Interpolate neutrino-nucleus pair annihilation kernals
!-----------------------------------------------------------------------

  fa0l             = fym * ( fdm * ( ftm * pa111  + ftp * pa121  )   &
&                  +         fdp * ( ftm * pa211  + ftp * pa221  ) ) &
&                  + fyp * ( fdm * ( ftm * pa112  + ftp * pa122  )   &
&                  +         fdp * ( ftm * pa212  + ftp * pa222  ) )

  fa0              = f10(fa0l)

  fp0l             = fym * ( fdm * ( ftm * pp111 + ftp * pp121 )   &
&                  +         fdp * ( ftm * pp211 + ftp * pp221 ) ) &
&                  + fyp * ( fdm * ( ftm * pp112 + ftp * pp122 )   &
&                  +         fdp * ( ftm * pp212 + ftp * pp222 ) )

  fp0              = f10(fp0l)

!-----------------------------------------------------------------------
!  Interpolate d(pair annihilation kernals)/d(rho)
!-----------------------------------------------------------------------

  fa0dl            = fdd * ( fym * ( ftm * ( -pa111  + pa211  )   &
&                  +                 ftp * ( -pa121  + pa221  ) ) &
&                  +         fyp * ( ftm * ( -pa112  + pa212  )   &
&                  +                 ftp * ( -pa122  + pa222  ) ) )

  fa0d             = ln_10 * fa0dl * fa0

  fp0dl            = fdd * ( fym * ( ftm * ( -pp111 + pp211 )   &
&                  +                 ftp * ( -pp121 + pp221 ) ) &
&                  +         fyp * ( ftm * ( -pp112 + pp212 )   &
&                  +                 ftp * ( -pp122 + pp222 ) ) )

  fp0d             = ln_10 * fp0dl * fp0

!-----------------------------------------------------------------------
!  Interpolate d(pair annihilation kernals)/d(t)
!-----------------------------------------------------------------------

  fa0tl            = ftt * ( fym * ( fdm * ( -pa111  + pa121  )   &
&                  +                 fdp * ( -pa211  + pa221  ) ) &
&                  +         fyp * ( fdm * ( -pa112  + pa122  )   &
&                  +                 fdp * ( -pa212  + pa222  ) ) )

  fa0t             = ln_10 * fa0tl * fa0

  fp0tl            = ftt * ( fym * ( fdm * ( -pp111 + pp121 )   &
&                  +                 fdp * ( -pp211 + pp221 ) ) &
&                  +         fyp * ( fdm * ( -pp112 + pp122 )   &
&                  +                 fdp * ( -pp212 + pp222 ) ) )

  fp0t             = ln_10 * fp0tl * fp0

!-----------------------------------------------------------------------
!  Interpolate d(pair annihilation kernals)/d(ye)
!-----------------------------------------------------------------------

  fa0yl            = fyy * ( fdm * ( ftm * ( -pa111  + pa112  )   &
&                  +                 ftp * ( -pa121  + pa122  ) ) &
&                  +         fdp * ( ftm * ( -pa211  + pa212  )   &
&                  +                 ftp * ( -pa221  + pa222  ) ) )

  fa0y             = ln_10 * fa0yl * fa0

  fp0yl            = fyy * ( fdm * ( ftm * ( -pp111 + pp112 )   &
&                  +                 ftp * ( -pp121 + pp122 ) ) &
&                  +         fdp * ( ftm * ( -pp211 + pp212 )   &
&                  +                 ftp * ( -pp221 + pp222 ) ) )

  fp0y             = ln_10 * fp0yl * fp0

!-----------------------------------------------------------------------
!  Construct neutrino-antineutrino pair annihilation functions
!-----------------------------------------------------------------------

  f0a (kp)         = fa0
  f0ad(kp)         = fa0d
  f0at(kp)         = fa0t
  f0ay(kp)         = fa0y
  f0p (kp)         = fp0
  f0pd(kp)         = fp0d
  f0pt(kp)         = fp0t
  f0py(kp)         = fp0y

END DO

RETURN
END SUBROUTINE pair_A_krnl

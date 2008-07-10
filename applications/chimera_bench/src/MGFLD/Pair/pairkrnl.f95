SUBROUTINE pairkrnl( n, j, ij_ray, ik_ray, k, rho, t, ye, f0a, f1a, f0ad, &
& f1ad, f0at, f1at, f0ay, f1ay, nez )
!-----------------------------------------------------------------------
!
!    File:         pairkrnl
!    Module:       pairkrnl
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
!  n              : neutrino type
!  j              : radial zone number
!  ij_ray         : j-index of a radial ray
!  ik_ray         : k-index of a radial ray
!  k              : neutrino energy grouip
!  nez            : neutrino energy array dimension
!  rho            : density (g/cm**3)
!  t              : temperature (K)
!  ye             : electron fraction
!
!    Output arguments:
!  f0a            : zero moment of the pair annihilation kernal
!  f01            : first moment of the pair annihilation kernal
!  f0ad           : derivative with respect to the density of the
!                    zero moment of the pair annihilation kernal
!  f1ad           : derivative with respect to the density of the
!                    first moment of the pair annihilation kernal
!  f0at           : derivative with respect to the temperature of the
!                    zero moment of the pair annihilation kernal
!  f1at           : derivative with respect to the temperature of the
!                    first moment of the pair annihilation kernal
!  f0ay           : derivative with respect to the electron fraction
!                    of the zero moment of the pair annihilation kernal
!  f1ay           : derivative with respect to the electron fraction
!                    of the first moment of the pair annihilation kernal
!
!    Variables that must be passed through common:
!
!  rhopairemn     : rho < rhopairemn: e-neutrino-antineutrino pair
!                    annihilation and production omitted
!  rhopairemx     : rho > rhopairemx: e-neutrino-antineutrino pair
!                    annihilation and production omitted
!  rhopairxmn     : rho < rhopairxmn: x-neutrino-antineutrino pair
!                    annihilation and production omitted
!  rhopairxmx     : rho > rhopairxmx: x-neutrino-antineutrino pair
!                    annihilation and production omitted
!   paira0i       : zero moment of the first neutrino-antineutrino
!                    pair annihilation kernal array
!   paira0ii      : zero moment of the second neutrino-antineutrino
!                    pair annihilation kernal array
!   paira1i       : first moment of the first neutrino-antineutrino
!                    pair annihilation kernal array
!   paira1ii      : first moment of the second neutrino-antineutrino
!                    pair annihilation kernal array
!  dgrid(idty(j,ij_ray,ik_ra))
!                 : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ra))
!                 : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ra))
!                 : number of table entries in ye between ye = 0.5
!                    and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)  : index for dgrid, tgrid, and ygrid for zone j
!   cv            : weak interaction coefficient
!   ca            : weak interaction coefficient
!  idrpp(j,ij_ray,ik_ray) : rho grid index for zone j
!  itrpp(j,ij_ray,ik_ray) : t grid index for zone j
!  iyrpp(j,ij_ray,ik_ray) : ye grid index for zone j
!
!    Subprograms called:
!      none
!
!    Include files:
!  array_module, kind_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, nu_energy_grid_module, pair_module,
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu, ij_ray_dim
USE numerical_module, ONLY : zero, half, one
USE physcnst_module, ONLY : cv, ca

USE edit_module, ONLY : nprint, nlog
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, idty
USE nu_energy_grid_module, ONLY : nnugp
USE pair_module, ONLY : idrpp, itrpp, iyrpp, paira0i, paira0ii
USE prb_cntl_module, ONLY : rhopairemn, rhopairemx, rhopairtmn, rhopairtmx

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

LOGICAL                            :: first = .true.
LOGICAL                            :: pair_off

INTEGER                            :: kp            ! antineutrino energy zone index

INTEGER                            :: id            ! EOS density index
INTEGER                            :: it            ! EOS temperature index
INTEGER                            :: iy            ! EOS electron fraction index

INTEGER                            :: i_ray         ! f(ij_ray,ik_ray)

INTEGER, PARAMETER                 :: ida = 1       ! lower cube table density index
INTEGER, PARAMETER                 :: idap1 = 2     ! upper cube table density index
INTEGER, PARAMETER                 :: ita = 1       ! lower cube table temperature index
INTEGER, PARAMETER                 :: itap1 = 2     ! upper cube table temperature index
INTEGER, PARAMETER                 :: iya = 1       ! lower cube table electron fraction index
INTEGER, PARAMETER                 :: iyap1 = 2     ! upper cube table electron fraction index

INTEGER, PARAMETER                 :: nnud=4        ! upper cube table electron fraction index

REAL(KIND=double), DIMENSION(nnud) :: cpair1       ! neutrino flavor coefficient
REAL(KIND=double), DIMENSION(nnud) :: cpair2       ! neutrino flavor coefficient

REAL(KIND=double)                  :: log_e         ! log10(e)
REAL(KIND=double)                  :: ln_10         ! ln(10)

REAL(KIND=double)                  :: fd            ! position of rho in grid
REAL(KIND=double)                  :: fdp           ! position of rho in grid wrt lower cube index
REAL(KIND=double)                  :: fdm           ! position of rho in grid wrt upper cube index
REAL(KIND=double)                  :: fdd           ! d(fd)/d(rho)
REAL(KIND=double)                  :: ft            ! position of t in grid
REAL(KIND=double)                  :: ftp           ! position of t in grid wrt lower cube index
REAL(KIND=double)                  :: ftm           ! position of t in grid wrt upper cube index
REAL(KIND=double)                  :: ftt           ! d(ft)/d(t)
REAL(KIND=double)                  :: fy            ! position of ye in grid
REAL(KIND=double)                  :: fyp           ! position of ye in grid wrt lower cube index
REAL(KIND=double)                  :: fym           ! position of ye in grid wrt upper cube index
REAL(KIND=double)                  :: fyy           ! d(fy)/d(ye)

REAL(KIND=double)                  :: pi111         ! scalar table entry for interpolation
REAL(KIND=double)                  :: pi211         ! scalar table entry for interpolation
REAL(KIND=double)                  :: pi121         ! scalar table entry for interpolation
REAL(KIND=double)                  :: pi112         ! scalar table entry for interpolation
REAL(KIND=double)                  :: pi221         ! scalar table entry for interpolation
REAL(KIND=double)                  :: pi212         ! scalar table entry for interpolation
REAL(KIND=double)                  :: pi122         ! scalar table entry for interpolation
REAL(KIND=double)                  :: pi222         ! scalar table entry for interpolation

REAL(KIND=double)                  :: pii111        ! scalar table entry for interpolation
REAL(KIND=double)                  :: pii211        ! scalar table entry for interpolation
REAL(KIND=double)                  :: pii121        ! scalar table entry for interpolation
REAL(KIND=double)                  :: pii112        ! scalar table entry for interpolation
REAL(KIND=double)                  :: pii221        ! scalar table entry for interpolation
REAL(KIND=double)                  :: pii212        ! scalar table entry for interpolation
REAL(KIND=double)                  :: pii122        ! scalar table entry for interpolation
REAL(KIND=double)                  :: pii222        ! scalar table entry for interpolation

REAL(KIND=double)                  :: j0il          ! interpolated log of zero-m pair funciton i
REAL(KIND=double)                  :: j0iil         ! interpolated log of zero-m pair funciton i
REAL(KIND=double)                  :: j0i           ! interpolated zero-m pair funciton i
REAL(KIND=double)                  :: j0ii          ! interpolated zero-m pair funciton i
REAL(KIND=double)                  :: j1i           ! interpolated first-m pair funciton i
REAL(KIND=double)                  :: j1ii          ! interpolated first-m pair funciton i

REAL(KIND=double)                  :: j0idl         ! interpolated log of d(j0i)/d(rho)
REAL(KIND=double)                  :: j0iidl        ! interpolated log of d(j0ii)/d(rho)
REAL(KIND=double)                  :: j0id          ! interpolated d(j0i)/d(rho)
REAL(KIND=double)                  :: j0iid         ! interpolated d(j0ii)/d(rho)
REAL(KIND=double)                  :: j1id          ! interpolated d(j1i)/d(rho)
REAL(KIND=double)                  :: j1iid         ! interpolated d(j1ii)/d(rho)

REAL(KIND=double)                  :: j0itl         ! interpolated log of d(j0i)/d(t)
REAL(KIND=double)                  :: j0iitl        ! interpolated log of d(j0ii)/d(t)
REAL(KIND=double)                  :: j0it          ! interpolated d(j0i)/d(t)
REAL(KIND=double)                  :: j0iit         ! interpolated d(jii)/d(t)
REAL(KIND=double)                  :: j1it          ! interpolated d(j1i)/d(t)
REAL(KIND=double)                  :: j1iit         ! interpolated d(j1ii)/d(t)

REAL(KIND=double)                  :: j0iyl         ! interpolated log of d(j0i)/d(ye)
REAL(KIND=double)                  :: j0iiyl        ! interpolated log of d(jii)/d(ye)
REAL(KIND=double)                  :: j0iy          ! interpolated d(j0i)/d(ye)
REAL(KIND=double)                  :: j0iiy         ! interpolated d(j0ii)/d(ye)
REAL(KIND=double)                  :: j1iy          ! interpolated d(j1i)/d(ye)
REAL(KIND=double)                  :: j1iiy         ! interpolated d(j1ii)/d(ye)

REAL(KIND=double)                  :: f10           ! 10** function
EXTERNAL f10

!-----------------------------------------------------------------------
!        Format
!-----------------------------------------------------------------------

 1001 FORMAT (' nnu=',i3,' exceeds nnud=',i3,' in subroutine pairkrnl')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) RETURN

!-----------------------------------------------------------------------
!  Compute standard model coefficients of the neutrino-antineutrino 
!   pair annihilation functions.
!-----------------------------------------------------------------------

IF (first) THEN
  first            = .false.

  IF ( nnu > nnud ) THEN
    WRITE (nprint,1001) nnu,nnud
    WRITE (nlog,1001) nnu,nnud
    STOP
  END IF

  cpair1(1)        = ( cv + ca )**2
  cpair2(1)        = ( cv - ca )**2
  IF ( nnu >= 2 ) THEN
    cpair1(2)      = ( cv - ca )**2
    cpair2(2)      = ( cv + ca )**2
  END IF
  IF ( nnu >= 3 ) THEN
    cpair1(3)      = ( cv + ca - 2.d0 )**2
    cpair2(3)      = ( cv - ca )**2
  END IF
  IF ( nnu >= 4 ) THEN
    cpair1(4)      = cpair2(3)
    cpair2(4)      = cpair1(3)
  END IF

  log_e            = DLOG10( DEXP( one ) )
  ln_10            = DLOG( 1.d+01 )

END IF

!-----------------------------------------------------------------------
!  Set pair annihilation functions to zero if rho outside specified
!   boundaries.
!-----------------------------------------------------------------------

pair_off           = .false.

IF ( n < 3  .and.  rho < rhopairemn ) pair_off = .true.
IF ( n < 3  .and.  rho > rhopairemx ) pair_off = .true.
IF ( n == 3  .and.  rho < rhopairtmn ) pair_off = .true.
IF ( n == 3  .and.  rho > rhopairtmx ) pair_off = .true.

IF ( pair_off ) THEN

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

END IF ! pair_off

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!  Compute independent variable grid indices and interpolation
!   coefficients.
!-----------------------------------------------------------------------

fd                 = dgrid(idty(j,ij_ray,ik_ray)) * DLOG10( rho )
id                 = idrpp(j,ij_ray,ik_ray)
fdp                = DMAX1( DMIN1( fd - DBLE( id ), one ), zero )
fdm                = one - fdp
fdd                = log_e * dgrid(idty(j,ij_ray,ik_ray))/( rho ) ! log_e = log(e)

ft                 = tgrid(idty(j,ij_ray,ik_ray)) * dlog10( t   )
it                 = itrpp(j,ij_ray,ik_ray)
ftp                = DMAX1( DMIN1( ft - DBLE( it ), one ), zero )
ftm                = one - ftp
ftt                = log_e * tgrid(idty(j,ij_ray,ik_ray))/( t   ) ! log_e = log(e)

fy                 = ygrid(idty(j,ij_ray,ik_ray)) * ( one - ye )
iy                 = iyrpp(j,ij_ray,ik_ray)
fyp                = DMAX1( DMIN1( fy - DBLE( iy ), one ), zero )
fym                = one - fyp
fyy                = - ygrid(idty(j,ij_ray,ik_ray))

DO kp = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Store neutrino-antineutrino pair annihilation functions in
!   temporary scalar variables for interpolation
!-----------------------------------------------------------------------

  pi111            = paira0i(j,k,kp,ida  ,ita  ,iya  ,i_ray)
  pi211            = paira0i(j,k,kp,idap1,ita  ,iya  ,i_ray)
  pi121            = paira0i(j,k,kp,ida  ,itap1,iya  ,i_ray)
  pi112            = paira0i(j,k,kp,ida  ,ita  ,iyap1,i_ray)
  pi221            = paira0i(j,k,kp,idap1,itap1,iya  ,i_ray)
  pi212            = paira0i(j,k,kp,idap1,ita  ,iyap1,i_ray)
  pi122            = paira0i(j,k,kp,ida  ,itap1,iyap1,i_ray)
  pi222            = paira0i(j,k,kp,idap1,itap1,iyap1,i_ray)

  pii111           = paira0ii(j,k,kp,ida  ,ita  ,iya  ,i_ray)
  pii211           = paira0ii(j,k,kp,idap1,ita  ,iya  ,i_ray)
  pii121           = paira0ii(j,k,kp,ida  ,itap1,iya  ,i_ray)
  pii112           = paira0ii(j,k,kp,ida  ,ita  ,iyap1,i_ray)
  pii221           = paira0ii(j,k,kp,idap1,itap1,iya  ,i_ray)
  pii212           = paira0ii(j,k,kp,idap1,ita  ,iyap1,i_ray)
  pii122           = paira0ii(j,k,kp,ida  ,itap1,iyap1,i_ray)
  pii222           = paira0ii(j,k,kp,idap1,itap1,iyap1,i_ray)

!-----------------------------------------------------------------------
!  Interpolate neutrino-antineutrino pair annihilation kernals
!-----------------------------------------------------------------------

  j0il             = fym * ( fdm * ( ftm * pi111  + ftp * pi121  )   &
&                  +         fdp * ( ftm * pi211  + ftp * pi221  ) ) &
&                  + fyp * ( fdm * ( ftm * pi112  + ftp * pi122  )   &
&                  +         fdp * ( ftm * pi212  + ftp * pi222  ) )

  j0i              = f10(j0il)
  j1i              = zero

  j0iil            = fym * ( fdm * ( ftm * pii111 + ftp * pii121 )   &
&                  +         fdp * ( ftm * pii211 + ftp * pii221 ) ) &
&                  + fyp * ( fdm * ( ftm * pii112 + ftp * pii122 )   &
&                  +         fdp * ( ftm * pii212 + ftp * pii222 ) )

  j0ii             = f10(j0iil)
  j1ii             = zero

!-----------------------------------------------------------------------
!  Interpolate d(pair annihilation kernals)/d(rho)
!-----------------------------------------------------------------------

  j0idl            = fdd * ( fym * ( ftm * ( -pi111  + pi211  )   &
&                  +                 ftp * ( -pi121  + pi221  ) ) &
&                  +         fyp * ( ftm * ( -pi112  + pi212  )   &
&                  +                 ftp * ( -pi122  + pi222  ) ) )

  j0id             = ln_10 * j0idl * j0i
  j1id             = zero

  j0iidl           = fdd * ( fym * ( ftm * ( -pii111 + pii211 )   &
&                  +                 ftp * ( -pii121 + pii221 ) ) &
&                  +         fyp * ( ftm * ( -pii112 + pii212 )   &
&                  +                 ftp * ( -pii122 + pii222 ) ) )

  j0iid            = ln_10 * j0iidl * j0ii
  j1iid            = zero

!-----------------------------------------------------------------------
!  Interpolate d(pair annihilation kernals)/d(t)
!-----------------------------------------------------------------------

  j0itl            = ftt * ( fym * ( fdm * ( -pi111  + pi121  )   &
&                  +                 fdp * ( -pi211  + pi221  ) ) &
&                  +         fyp * ( fdm * ( -pi112  + pi122  )   &
&                  +                 fdp * ( -pi212  + pi222  ) ) )

  j0it             = ln_10 * j0itl * j0i
  j1it             = zero

j0iitl             = ftt * ( fym * ( fdm * ( -pii111 + pii121 )   &
&                  +                 fdp * ( -pii211 + pii221 ) ) &
&                  +         fyp * ( fdm * ( -pii112 + pii122 )   &
&                  +                 fdp * ( -pii212 + pii222 ) ) )

  j0iit            = ln_10 * j0iitl * j0ii
  j1iit            = zero

!-----------------------------------------------------------------------
!  Interpolate d(pair annihilation kernals)/d(ye)
!-----------------------------------------------------------------------

  j0iyl            = fyy * ( fdm * ( ftm * ( -pi111  + pi112  )   &
&                  +                 ftp * ( -pi121  + pi122  ) ) &
&                  +         fdp * ( ftm * ( -pi211  + pi212  )   &
&                  +                 ftp * ( -pi221  + pi222  ) ) )

  j0iy             = ln_10 * j0iyl * j0i
  j1iy             = zero

  j0iiyl           = fyy * ( fdm * ( ftm * ( -pii111 + pii112 )   &
&                  +                 ftp * ( -pii121 + pii122 ) ) &
&                  +         fdp * ( ftm * ( -pii211 + pii212 )   &
&                  +                 ftp * ( -pii221 + pii222 ) ) )

  j0iiy            = ln_10 * j0iiyl * j0ii
  j1iiy            = zero

!-----------------------------------------------------------------------
!  Construct neutrino-antineutrino pair annihilation functions
!-----------------------------------------------------------------------

  f0a (kp)         = cpair1(n) * j0i  + cpair2(n) * j0ii
  f0ad(kp)         = cpair1(n) * j0id + cpair2(n) * j0iid
  f0at(kp)         = cpair1(n) * j0it + cpair2(n) * j0iit
  f0ay(kp)         = cpair1(n) * j0iy + cpair2(n) * j0iiy
  f1a (kp)         = cpair1(n) * j1i  + cpair2(n) * j1ii
  f1ad(kp)         = cpair1(n) * j1id + cpair2(n) * j1iid
  f1at(kp)         = cpair1(n) * j1it + cpair2(n) * j1iit
  f1ay(kp)         = cpair1(n) * j1iy + cpair2(n) * j1iiy

END DO

RETURN
END SUBROUTINE pairkrnl

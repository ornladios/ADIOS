SUBROUTINE pair_A_ed( n, j, ij_ray, ik_ray, k, r, rho, t, ye, a0w, b0w, &
& c0w, a1w, b1w, c1w, rmdnts0, rmdnts1, rmdnts, artpe, brtpe, crtpe,    &
& artae, brtae, crtae )
!-----------------------------------------------------------------------
!
!    File:         pair_A_ed
!    Module:       pair_A_ed
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/07/03
!
!    Purpose:c
!      To compute quantities which are needed for editing neutrino-antineutrino
!       pair annihilation. In particular the quantities rmdnts0, rmdnts1, and
!       rmdnts are computed, which are the inverse mean free paths for
!       neutrino-antineutrino pair annihilation. The quantities artpe, brtpe,
!       crtpe, artae, brtae, are crtae from which the pair production and
!       annihilation rates, and the mean neutrino energy annihilated and
!       created can be computed. The quantities a0w, b0w, c0w, a1w, b1w,
!       and c1w are computed from which zero and first Legendre coefficients
!       for the n-type neutrino pair annihilation and production functions
!       These are included in the multi-group diffusion equations, which have
!       the form
!
!            pis1 = -yt*( dpsi0/dr - a1w*psi0 - c1w)
!
!            dpsi0/dt/! + d(r2*psi1)/dr/3./r2 = xw + yw*psi0 + zw*dpsi0/dr
!
!         where
!
!            1/yt = zltot = jw + yaw - b1w
!            xw   = jw + c0w + b0w*c1w*yt
!            yw   = -jw - yaw + a0w + b0w*a1w*yt
!            zw   = -b0w*yt
!
!         The quantities a0w, b0w, c0w, a1w, b1w, c1w are given by
!
!        a0w  = -K   Int{ w2'dw'[ phi0p(w,w')( 1 - psi0bar(w') ) + phi0a(w,w')psi0bar(w') ] }
!        b0w  =  K/3 Int{ w2'dw'[ phi1p(w,w') - phi1a(w,w') ]psi1bar(w') }
!        c0w  =  K   Int{ w2'dw'[ phi0p(w,w')( 1 - psi0bar(w') ) ] }
!        a1w  =  K   Int{ w2'dw'[ phi1p(w,w') - phi1a(w,w') ]psi1bar(w') }
!        b1w  = -K   Int{ w2'dw'[ phi0p(w,w')( 1 - psi0bar(w') ) + phi0a(w,w')psi0bar(w') ] }
!        c1w  = -K   Int{ w2'dw'[ phi1p(w,w')psi1bar(w') ] }
!
!         where
!
!            K    = 2*pi/c * 1/(hc)**3
!
!        The Legendre moments of the neutrino pair annihilation and production
!         functions are given by
!
!          phiLa  = ( g2/pi ) * cpair1(n)*jaLi(w,w') + cpair2(n)*jaLii(w,w')
!
!        where
!
!      jaLi(w,w') = Int{ de*[1 - Fe(e)]*[1 - F(w' + w - e)] * raLi(e,w,w') }
!
!        and jaLii(w,w') is defined likewise.
!
!          phiLp  = ( g2/pi ) * cpair1(n)*jpLi(w,w') + cpair2(n)*jpLii(w,w')
!
!        where
!
!      jpLi(w,w') = Int{ de*Fe(e)*F(w' + w - e) * rpLi(e,w,w') }
!
!        and jpLii(w,w') is defined likewise.
!
!        The integrations are performed in subroutine paircal.
!
!      Finally, the quantities artpe, brtpe, crtpe, artae, brtae, and crtae
!       are computed here. These are used for computing the production rate
!       of n-type neutrinos by electron- positron pair annihilation and of the
!       annihilation rate of n-type neutrinos by the inverse process.
!
!      The zero moment of the production rate per state, dpsi0dt, of neutrinos
!       of type n is given by
!
!           dppsi0dt = dppsi0/dtnph
!
!       where
!
!            dppis0  = c*[ artpe*psi0 + brtpe*psi1 + crtpe ]
!
!       and where
!
!      artpe  = -K   Int{ w2'dw'[ phi0p(w,w')( 1 - psi0bar(w') ) ] }
!      brtpe  =  K/3 Int{ w2'dw'[ phi1p(w,w')psi1bar(w') ] }
!      crtpe  =  K   Int{ w2'dw'[ phi0p(w,w')( 1 - psi0bar(w') ) ] }
!
!      The zero moment of the annihilation rate per state, dapsi0dt,
!       of neutrinos of type n is given by
!
!           dapsi0dt = dapsi0/dtnph
!
!       where
!
!            dapis0  = !*[ artae*psi0 + brtae*psi1 + crtae ]
!
!       and where
!
!      artae  =  K   Int{ w2'dw'[ phi0a(w,w')psi0bar(w') ] }
!      brtae  = aK/3 Int{ w2'dw'[ phi1a(w,w')psi1bar(w') ] }
!      crtae  =  0.0
!
!       (Note that dapis0 is computed to be positive so that
!        a positive rate results.)
!
!    Input arguments:
!
!  n             : neutrino type
!  j             : radial zone number
!  ij_ray        : j-index of a radial ray
!  ik_ray        : k-index of a radial ray
!  k             : energy zone number
!  r             : radius (cm)
!  rho           : matter density (g/cm**3)
!  t             : matter temperature (K)
!  ye            : electron fraction
!
!    Output arguments:
!
!  a0w           : coefficient for computing net pair annihilation rates
!  b0w           : coefficient for computing net pair annihilation rates
!  c0w           : coefficient for computing net pair annihilation rates
!  a1w           : coefficient for computing net pair annihilation rates
!  b1w           : coefficient for computing net pair annihilation rates
!  c1w           : coefficient for computing net pair annihilation rates
!
!  artpe, brtpe, crtpe, artae, brtae, and crtae
!                : coefficients from which the rate of neutrino-antineutrino
!                   pair production rate can be computed. These coefficients
!                   are used in subroutine editng.
!
!    Input arguments (common):
!
!  nnugp(n)      : number of energy zones for neutrinos of type n
!  unu(k)        : energy of energy zone k (MeV)
!  dunu(k)       : energy width of energy zone k (MeV)
!  psi0(j,k,n)   : zeroth moment of of the neutrino occupation probability
!                   for neutrinos of type n,  energy zone k, radial zone j
!  psi1(j,k,n)   : first moment of of the neutrino occupation probability for
!                   neutrinos of type n, energy zone k, radial zone j
!  r_sphere(k,n) : radius of the neutrinosphere for neutrinos of type n and
!                   energy k (cm).
!
!    Subprograms called:
!  pairkrnl      : interpolates the pair annihilation functions
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, nu_dist_module, nu_energy_grid_module, pair_A_module,
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY : nez, nnu
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : kmev

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : unu, dunu, psi0, psi1, r_sphere
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE pair_A_module
USE prb_cntl_module, ONLY : ipairA, rhopairAemn, rhopairAemx, rhopairAtmn, rhopairAtmx

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

REAL(KIND=double), INTENT(in)     :: r             ! radius of the neutrino-sphere (cm)
REAL(KIND=double), INTENT(in)     :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)     :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)     :: ye            ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)    :: rmdnts        ! neutrino inverse mean free path
REAL(KIND=double), INTENT(out)    :: rmdnts0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double), INTENT(out)    :: rmdnts1       ! first moment of the neutrino inverse mean free path

REAL(KIND=double), INTENT(out)    :: a0w           ! coefficient for computing change in psi0
REAL(KIND=double), INTENT(out)    :: b0w           ! coefficient for computing change in psi0
REAL(KIND=double), INTENT(out)    :: c0w           ! coefficient for computing change in psi0
REAL(KIND=double), INTENT(out)    :: a1w           ! coefficient for computing change in psi1
REAL(KIND=double), INTENT(out)    :: b1w           ! coefficient for computing change in psi1
REAL(KIND=double), INTENT(out)    :: c1w           ! coefficient for computing change in psi1

REAL(KIND=double), INTENT(out)    :: artpe         ! coefficient for computing pair production rate
REAL(KIND=double), INTENT(out)    :: brtpe         ! coefficient for computing pair production rate
REAL(KIND=double), INTENT(out)    :: crtpe         ! coefficient for computing pair production rate
REAL(KIND=double), INTENT(out)    :: artae         ! coefficient for computing pair annihilation rate
REAL(KIND=double), INTENT(out)    :: brtae         ! coefficient for computing pair annihilation rate
REAL(KIND=double), INTENT(out)    :: crtae         ! coefficient for computing pair annihilation rate

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                :: var_name

LOGICAL                           :: pair_off
LOGICAL                           :: first = .true.

INTEGER                           :: kp            ! outcoming neutrino energy zone index

INTEGER, DIMENSION(4)             :: nanti         ! neutrino pair identifier
INTEGER                           :: na            ! antineutrino flavor index

INTEGER                           :: istat         ! allocation status

REAL(KIND=double)                 :: enu           ! neutrino energy/kt
REAL(KIND=double)                 :: enubar        ! antineutrino energy/kt

REAL(KIND=double)                 :: x1            ! psi1
REAL(KIND=double)                 :: x0a           ! psi1
REAL(KIND=double)                 :: x1a           ! psi1

REAL(KIND=double)                 :: c1            ! integrated scattering function

REAL(KIND=double)                 :: art           ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: brt           ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: drt           ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: ert           ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: zlthp         ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: asthp         ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: gthp          ! coefficient for computing the mean pair-annihilation rate

REAL(KIND=double), PARAMETER      :: rsinmin = .99d+00 ! maximum value of r_sphere/r
REAL(KIND=double)                 :: psib10        ! psi1/psi0 for antineutrinos
REAL(KIND=double)                 :: psi10         ! psi1/psi0 for neutrinos
REAL(KIND=double)                 :: ac1           ! 1 - psi10 * psib10/6
REAL(KIND=double)                 :: rsphdr        ! R(neutrino-sphere)/radius
REAL(KIND=double)                 :: rsphbdr       ! R(antineutrino-sphere)/radius
REAL(KIND=double)                 :: z             ! sqrt(1 - rsphdr**2)
REAL(KIND=double)                 :: zb            ! sqrt(1 - rsphbdr**2)
REAL(KIND=double)                 :: zav           ! ( 1 + z  )/2
REAL(KIND=double)                 :: zbav          ! ( 1 + zb  )/2
REAL(KIND=double)                 :: z2av          ! ( 1 + z + z**2 )/3
REAL(KIND=double)                 :: zb2av         ! ( 1 + zb + zb**2 )/3
REAL(KIND=double)                 :: ac2           ! cuttoff for pair-annihilation rate due to radial streaming
REAL(KIND=double)                 :: ac            ! minimum of ac1 and ac2

REAL(KIND=double)                 :: fexp          ! exponential function


REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: w2dw         ! w2 * neutrino energy bin width

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0a          ! zero moment of nuclear the pair-annihilation function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ad         ! d(f0a)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0at         ! d(f0a)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ay         ! d(f0a)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1a          ! first moment of the nuclear pair-annihilation function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ad         ! d(f1a)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1at         ! d(f1a)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ay         ! d(f1a)/d(ye)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0p          ! zero moment of the nuclear pair-production function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0pd         ! d(f0a)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0pt         ! d(f0a)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0py         ! d(f0a)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1p          ! first moment of the nuclear pair-production function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1pd         ! d(f1a)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1pt         ! d(f1a)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1py         ! d(f1a)/d(ye)

EXTERNAL fexp

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in paired')
 2001 FORMAT (' Deallocation problem for array ',a10,' in paired')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!  No pair annihilation and production if ipairA = 0
!-----------------------------------------------------------------------


pair_off             = .false.
IF ( ipairA == 0 ) pair_off = .true.
IF ( nnugp(n) == 0 ) pair_off = .true.
IF ( n <  3  .and.  rho < rhopairAemn ) pair_off = .true.
IF ( n <  3  .and.  rho > rhopairAemx ) pair_off = .true.
IF ( n == 3  .and.  rho < rhopairAtmn ) pair_off = .true.
IF ( n == 3  .and.  rho > rhopairAtmx ) pair_off = .true.
IF ( pair_off ) THEN

  a0w                = zero
  b0w                = zero
  c0w                = zero
  a1w                = zero
  b1w                = zero
  c1w                = zero
  rmdnts0            = zero
  rmdnts1            = zero
  rmdnts             = zero
  artpe              = zero
  brtpe              = zero
  crtpe              = zero
  artae              = zero
  brtae              = zero
  crtae              = zero

  RETURN

END IF ! pair_off

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (w2dw(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w2dw      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (f0a(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0a       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0ad(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ad      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0at(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0at      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0ay(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ay      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1a(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1a       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1ad(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ad      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1at(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1at      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1ay(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ay      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (f0p(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0p       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0pd(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0pd      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0pt(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0pt      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f0py(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0py      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1p(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1p       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1pd(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1pd      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1pt(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1pt      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (f1py(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1py      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  nanti(1)         = 2
  nanti(2)         = 1
  nanti(3)         = 4
  nanti(4)         = 3
END IF ! first

w2dw               = zero

f0a                = zero
f0ad               = zero
f0at               = zero
f0ay               = zero
f1a                = zero
f1ad               = zero
f1at               = zero
f1ay               = zero

f0p                = zero
f0pd               = zero
f0pt               = zero
f0py               = zero
f1p                = zero
f1pd               = zero
f1pt               = zero
f1py               = zero

!-----------------------------------------------------------------------
!  Compute variables related to n-type neutrinos.
!
!  tmev       : t (MeV)
!  enu        : neutrino energy (MeV)
!  enubar     : antineutrino energy (MeV)
!  na         : neutrino index corresponding to the antiparticle of the neutrino corresponding to index n
!-----------------------------------------------------------------------

enu                = unu(j,k)
na                 = nanti(n)

!-----------------------------------------------------------------------
!       ||||| Computation of n-neutrino pair annihilation |||||
!       |||||          and production functions           |||||
!-----------------------------------------------------------------------

zlthp              = zero
asthp              = zero
c1                 = zero
gthp               = zero
art                = zero
brt                = zero
drt                = zero
ert                = zero

DO kp = 1,nnugpmx
  w2dw(kp)         = unu(j,kp)**2 * dunu(j,kp)
END DO

!-----------------------------------------------------------------------
!  Get pair annihilation kernels
!-----------------------------------------------------------------------

CALL pair_A_krnl( n, j, ij_ray, ik_ray, k, rho, t, ye, f0a, f0p, f0ad, &
& f0pd, f0at, f0pt, f0ay, f0py, nez )

x1                 = half * ( psi1(j-1,k,n ) + psi1(j,k,n ) )

!-----------------------------------------------------------------------
!                   ||||| Begin loop over kp |||||
!-----------------------------------------------------------------------

DO kp = 1,nnugp(n)

  enubar           = unu(j,kp)
  x0a              = psi0(j,kp,na)
  x1a              = half * ( psi1(j-1,kp,na) + psi1(j,kp,na) )

!-----------------------------------------------------------------------
!  Compute angular correction factors to  pair annihilation kernals
!-----------------------------------------------------------------------

  psib10           = DMIN1( DABS( x1a)/( x0a + epsilon ), one )
  psi10            = DMIN1( DABS( x1 )/( psi0(j,k,n ) + epsilon ), one )
  ac1              = one - psi10 * psib10/6.d0
  rsphdr           = DMIN1( r_sphere(k ,n ,ij_ray,ik_ray)/r, rsinmin )
  z                = DSQRT( one - rsphdr**2  )
  rsphbdr          = DMIN1( r_sphere(kp,na,ij_ray,ik_ray)/r, rsinmin )
  zb               = DSQRT( one - rsphbdr**2 )
  zav              = ( one + z  )/2.d0
  zbav             = ( one + zb )/2.d0
  z2av             = ( one + z  + z**2  )/3.d0
  zb2av            = ( one + zb + zb**2 )/3.d0
  ac2              = 0.75d0 * ( one - 2.d0 * zav * zbav + z2av * zb2av + 0.5d0 * ( one - z2av ) * ( one - zb2av ) )
  ac               = DMAX1( ac1, ac2, zero )

!-----------------------------------------------------------------------
!  Apply angular correction to pair annihilation kernals
!-----------------------------------------------------------------------

  f0a(kp)          = f0a(kp) * ac
  f1a(kp)          = f1a(kp) * ac

!-----------------------------------------------------------------------
!  Assemble pair annihilation and production functions
!-----------------------------------------------------------------------

  zlthp            = zlthp + ( f0a(kp) * x0a + f0p(kp) * ( one - x0a ) )   * w2dw(kp)
  asthp            = asthp + x1a * ( f1a(kp) - f1p(kp) )   * w2dw(kp)
  c1               = c1    + f0p(kp) * ( one - x0a )       * w2dw(kp)
  gthp             = gthp  - x1a * f1p(kp)                 * w2dw(kp)
  art              = art   + f0p(kp) * ( one - x0a )       * w2dw(kp)
  brt              = brt   + f1p(kp) * x1a                 * w2dw(kp)
  drt              = drt   + f0a(kp) * x0a                 * w2dw(kp)
  ert              = ert   + f1a(kp) * x1a                 * w2dw(kp)

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

END DO

rmdnts0            =  zlthp
rmdnts             =  zlthp
a0w                = -zlthp
b0w                = -asthp/3.d0
c0w                =  c1
a1w                = -asthp
b1w                = -zlthp
c1w                =  gthp
artpe              = -art
brtpe              =  brt/3.d0
crtpe              =  art
artae              =  drt
brtae              =  ert/3.d0
crtae              =  zero

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (w2dw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w2dw      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (f0a, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0a       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0ad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ad      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0at, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0at      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0ay, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0ay      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1a, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1a       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1ad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ad      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1at, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1at      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1ay, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ay      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (f0p, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0p       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0pd, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0pd      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0pt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0pt      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f0py, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f0py      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1p, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1p       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1pd, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1pd      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1pt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1pt      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1py, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1py      '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE pair_A_ed

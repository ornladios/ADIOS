SUBROUTINE bremed( n, j, ij_ray, ik_ray, k, r, rho, t, ye, a0w, b0w, c0w, &
& a1w, b1w, c1w, rmdnts0, rmdnts1, rmdnts, artpe, brtpe, crtpe, artae, &
& brtae, crtae )
!-----------------------------------------------------------------------
!
!    File:         bremed
!    Module:       bremed
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/08/02
!
!    Purpose:
!      To compute quantities which are needed for editing neutrino-antineutrino
!       bremsstrahlung pair annihilation. In particular the quantity rmdnts is
!       computed, which is the inverse mean free path for neutrino-antineutrino
!       bremsstrahlung pair annihilation. The quantities a0w, b0w, c0w, a1w, b1w,
!       and c1w are computed from which zero and first Legendre coefficients for
!       the n-type neutrino bremsstrahlung pair annihilation and production
!       functions. These are included in the multi-group diffusion equations,
!       which have the form
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
!            K    = 2*pi/! * 1/(hc)**3
!
!        The Legendre moments of the neutrino bremsstrahlung pair annihilation
!         and production functions are given by
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
!      Finally, the quantities artpe, brtpe, crtpe, artae, brtae, and crtae are
!       computed here. These are used for computing the production rate of n-type
!       neutrinos by electron- positron pair annihilation and of the annihilation
!       rate of n-type neutrinos by the inverse process.
!
!      The zero moment of the production rate per state, dpsi0dt, of neutrinos
!       of type n is given by
!
!           dppsi0dt = dppsi0/dtnph
!
!       where
!
!            dppis0  = !*[ artpe*psi0 + brtpe*psi1 + crtpe ]
!
!       and where
!
!      artpe  = -K   Int{ w2'dw'[ phi0p(w,w')( 1 - psi0bar(w') ) ] }
!      brtpe  =  K/3 Int{ w2'dw'[ phi1p(w,w')psi1bar(w') ] }
!      crtpe  =  K   Int{ w2'dw'[ phi0p(w,w')( 1 - psi0bar(w') ) ] }
!
!      The zero moment of the annihilation rate per state, dapsi0dt, of neutrinos of type n is
!       given by
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
!       (Note that dapis0 is computed to be positive so that a positive rate results.)
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
!                : coefficients from which the rate of neutrino-
!          antineutrino pair production rate can be computed. These
!          coefficients are used in subroutine editng.
!
!    Input arguments (common):
!
!  nnugp(n)      : number of energy zones for neutrinos of type n
!  unu(k)        : energy of energy zone k (MeV)
!  dunu(k)       : energy width of energy zone k (MeV)
!  psi0(j,k,n)   : zeroth moment of of the neutrino occupation probability for
!                   neutrinos of type n, energy zone k, radial zone j
!  psi1(j,k,n)   : first moment of of the neutrino occupation probability for
!                   neutrinos of type n, energy zone k, radial zone j
!  r_sphere(k,n) : radius of the neutrinosphere for neutrinos of type n and energy k (cm).
!
!    Subprograms called:
!      bremrkrnl
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, nu_dist_module, nu_energy_grid_module, prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nez, nnu
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : kmev

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : unu, dunu, psi0, psi1, r_sphere
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : ibrem, rhobrememn, rhobrememx, rhobremtmn, rhobremtmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: j             ! radial zone index
INTEGER, INTENT(in)               :: k             ! incoming neutrino energy zone index
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)     :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)     :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)     :: ye            ! electron fraction
REAL(KIND=double), INTENT(in)     :: r             ! radius (cm)

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

REAL(KIND=double), INTENT(out)    :: artpe         ! coefficient for computing B pair production rate
REAL(KIND=double), INTENT(out)    :: brtpe         ! coefficient for computing B pair production rate
REAL(KIND=double), INTENT(out)    :: crtpe         ! coefficient for computing B pair production rate
REAL(KIND=double), INTENT(out)    :: artae         ! coefficient for computing B pair annihilation rate
REAL(KIND=double), INTENT(out)    :: brtae         ! coefficient for computing B pair annihilation rate
REAL(KIND=double), INTENT(out)    :: crtae         ! coefficient for computing B pair annihilation rate

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                :: var_name

LOGICAL                           :: brem_off
LOGICAL                           :: first = .true.

INTEGER                           :: kb            ! antineutrino energy zone index

INTEGER, DIMENSION(4)             :: nanti
INTEGER                           :: na            ! antineutrino flavor index

INTEGER                           :: istat         ! allocation status

REAL(KIND=double)                 :: f0p           ! pair production kernel
REAL(KIND=double)                 :: f1p           ! pair production kernel

REAL(KIND=double)                 :: art           ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: brt           ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: drt           ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: ert           ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: zlthp         ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: asthp         ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: gthp          ! coefficient for computing the mean pair-annihilation rate
REAL(KIND=double)                 :: coefpp        ! coefficient for computing prodction from annihilation rates

REAL(KIND=double)                 :: enu           ! neutrino energy/kt
REAL(KIND=double)                 :: enubar        ! antineutrino energy/kt
REAL(KIND=double)                 :: tmev          ! temperature (MeV)

REAL(KIND=double)                 :: x1            ! psi1
REAL(KIND=double)                 :: x0a           ! psi1
REAL(KIND=double)                 :: x1a           ! psi1

REAL(KIND=double)                 :: c1            ! integrated B pair production function

REAL(KIND=double), PARAMETER      :: rsinmin = .99d+00 ! maximum value of r_sphere/r
REAL(KIND=double)                 :: ac1           ! 1 - psi10 * psib10/6
REAL(KIND=double)                 :: ac2           ! cuttoff for pair-annihilation rate due to radial streaming
REAL(KIND=double)                 :: ac            ! minimum of ac1 and ac2
REAL(KIND=double)                 :: z             ! sqrt(1 - rsphdr**2)
REAL(KIND=double)                 :: zb            ! sqrt(1 - rsphbdr**2)
REAL(KIND=double)                 :: zav           ! ( 1 + z  )/2
REAL(KIND=double)                 :: zbav          ! ( 1 + zb  )/2
REAL(KIND=double)                 :: z2av          ! ( 1 + z + z**2 )/3
REAL(KIND=double)                 :: zb2av         ! ( 1 + zb + zb**2 )/3
REAL(KIND=double)                 :: psib10        ! psi1/psi0 for antineutrinos
REAL(KIND=double)                 :: psi10         ! psi1/psi0 for neutrinos
REAL(KIND=double)                 :: rsphdr        ! R(neutrino-sphere)/radius
REAL(KIND=double)                 :: rsphbdr       ! R(antineutrino-sphere)/radius

REAL(KIND=double)                 :: fexp          ! exponential function

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: w2dw         ! w2 * neutrino energy bin width

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0a          ! zero moment of the pair-annihilation function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ad         ! d(f0a)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0at         ! d(f0a)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f0ay         ! d(f0a)/d(ye)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1a          ! first moment of the pair-annihilation function
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ad         ! d(f1a)/d(rho)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1at         ! d(f1a)/d(t)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: f1ay         ! d(f1a)/d(ye)

EXTERNAL fexp

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for bremred arrays')
 2001 FORMAT (' Deallocation problem for bremred arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set bremsstrahlung pair annihilation functions to zero if
!     nnugp(n) =  0
!   or
!     ibrem = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0  .or.  ibrem == 0 ) THEN

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

END IF ! ibrem = 0

!-----------------------------------------------------------------------
!  Set bremsstrahlung pair annihilation functions to zero if rho
!   outside specified boundaries.
!-----------------------------------------------------------------------

brem_off             = .false.

IF ( n <  3  .and.  rho < rhobrememn ) brem_off = .true.
IF ( n <  3  .and.  rho > rhobrememx ) brem_off = .true.
IF ( n >= 3  .and.  rho < rhobremtmn ) brem_off = .true.
IF ( n >= 3  .and.  rho > rhobremtmx ) brem_off = .true.

IF ( brem_off ) THEN
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
END IF ! brem_off

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

!-----------------------------------------------------------------------
!  Compute variables related to n-type neutrinos.
!
!  tmev       : t (MeV)
!  enu        : neutrino energy (MeV)
!  enubar     : antineutrino energy (MeV)
!  na         : neutrino index corresponding to the antiparticle of
!                the neutrino corresponding to index n
!-----------------------------------------------------------------------

enu                = unu(j,k)
tmev               = t * kmev
na                 = nanti(n)

!-----------------------------------------------------------------------
!       ||||| Computation of n-neutrino bremsstrahlung pair |||||
!       |||||    annihilation and production functions      |||||
!-----------------------------------------------------------------------

zlthp              = zero
asthp              = zero
c1                 = zero
gthp               = zero
art                = zero
brt                = zero
drt                = zero
ert                = zero

DO kb = 1,nnugpmx
  w2dw(kb)         = unu(j,kb)**2 * dunu(j,kb)
END DO

!-----------------------------------------------------------------------
!  Get pair annihilation kernels
!-----------------------------------------------------------------------

CALL bremkrnl( n, j, ij_ray, ik_ray, k, rho, t, ye, f0a, f1a, f0ad, f1ad, &
& f0at, f1at, f0ay, f1ay, nez )

DO kb = 1,nnugp(n)

  enubar           = unu(j,kb)
  x0a              = psi0(j,kb,na)
  x1a              = half * ( psi1(j-1,kb,na) + psi1(j,kb,na) )
  x1               = half * ( psi1(j-1,kb,n ) + psi1(j,kb,n ) )

  coefpp           = fexp( -( enu + enubar )/tmev )

!-----------------------------------------------------------------------
!  Compute pair production kernals using detailed balance
!-----------------------------------------------------------------------

  f0p              = coefpp * f0a(kb)
  f1p              = coefpp * f1a(kb)

!-----------------------------------------------------------------------
!  Compute angular correction factors to bremsstrahlung pair
!   annihilation kernals
!-----------------------------------------------------------------------

  psib10           = DABS( x1a)/( x0a + epsilon )
  psi10            = DABS( x1 )/( psi0(j,k,n ) + epsilon )
  ac1              = one - psi10 * psib10/6.d0
  rsphdr           = DMIN1( r_sphere(k,n ,ij_ray,ik_ray)/r, rsinmin )
  z                = DSQRT( one - rsphdr**2  )
  rsphbdr          = DMIN1( r_sphere(kb,na,ij_ray,ik_ray)/r, rsinmin )
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

  f0a(kb)          = f0a(kb) * ac
  f1a(kb)          = f1a(kb) * ac

!-----------------------------------------------------------------------
!  Assemble bremsstrahlung pair annihilation and production functions
!-----------------------------------------------------------------------

  zlthp            = zlthp + ( f0a(kb) * x0a + f0p * ( one - x0a ) )   * w2dw(kb)
  asthp            = asthp + x1a * ( f1a(kb) - f1p )                   * w2dw(kb)
  c1               = c1    + f0p * ( one - x0a )                       * w2dw(kb)
  gthp             = gthp  - x1a * f1p                                 * w2dw(kb)
  art              = art   + f0p * ( one - x0a )                       * w2dw(kb)
  brt              = brt   + f1p * x1a                                 * w2dw(kb)
  drt              = drt   + f0a(kb) * x0a                             * w2dw(kb)
  ert              = ert   + f1a(kb) * x1a                             * w2dw(kb)

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

RETURN
END

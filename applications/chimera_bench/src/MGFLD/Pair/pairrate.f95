SUBROUTINE pairrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, r, &
& nx )
!-----------------------------------------------------------------------
!
!    File:         pairrate
!    Module:       pairrate
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/07/03
!
!    Purpose:
!      The zero and first legendre coefs for the n-type neutrino pair
!       annihilation and  production functions are computed here.
!      These are included in the multi-group diffusion equations,
!       which have the form
!
!            pis1 = -yt*( dpsi0/dr - a1w*psi0 - c1w)
!
!            dpsi0/dt/! + d(r2*psi1)/dr/3./r2 = xw + yw*psi0 + zw*dpsi0/dr
!
!       where
!
!            1/yt = zltot = jw + yaw - b1w
!            xw   = jw + c0w + b0w*c1w*yt
!            yw   = -jw - yaw + a0w + b0w*a1w*yt
!            zw   = -b0w*yt
!
!      Neutrino pair annihilation and production contributes to the
!       terms a0w, a1w, b0w, b1w, c0w, c1w as follows:
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
!      The Legendre moments of the neutrino pair annihilation and production
!       functions are given by
!
!          phiLa  = ( g2/pi ) * cpair1(n)*jaLi(w,w') + cpair2(n)*jaLii(w,w')
!
!       where
!
!      jaLi(w,w') = Int{ de*[1 - Fe(e)]*[1 - F(w' + w - e)] * raLi(e,w,w') }
!
!       and jaLii(w,w') is defined likewise.
!
!          phiLp  = ( g2/pi ) * cpair1(n)*jpLi(w,w') + cpair2(n)*jpLii(w,w')
!
!       where
!
!      jpLi(w,w') = Int{ de*Fe(e)*F(w' + w - e) * rpLi(e,w,w') }
!
!       and jpLii(w,w') is defined likewise.
!
!          phiLp  = ( g2/pi ) * cpair1(n)*jpLi(w,w') + cpair2(n)*jpLii(w,w')
!
!       where
!
!      jpLi(w,w') = Int{ de*Fe(e)*F(w' + w - e) * rpLi(e,w,w') }
!
!       and jpLii(w,w') is defined likewise.
!
!      The integrations are performed in subroutine paircal.
!
!    Input arguments:
!
!  n             : neutrino type
!  jr_min        : inner radial zone number
!  jr_max        : outer radial zone number
!  ij_ray        : j-index of a radial ray
!  ik_ray        : k-index of a radial ray
!
!    Input arguments (common):
!
!  r(j)          : radius (cm)
!  rho(j)        : matter density (g/cm**3)
!  t(j)          : matter temperature (K)
!  ye(j)         : electron fraction
!  ipair         : 0, neutrino pair annihilation and production turned off;
!                   pair annihilation subroutines are bypassed; pair annihilation
!                   function arrays, if used, must be zeroed elsewhere
!                  1, neutrino pair annihilation and production included
!  rhopairemn    : rho < rhopairemn: e-neutrino-antineutrino pair
!                   annihilation and production omitted
!  rhopairemx    : rho > rhopairemx: e-neutrino-antineutrino pair
!                   annihilation and production omitted
!  rhopairxmn    : rho < rhopairxmn: x-neutrino-antineutrino pair
!                   annihilation and production omitted
!  rhopairxmx    : rho > rhopairxmx: x-neutrino-antineutrino pair
!                   annihilation and production omitted
!  nnugp(n)      : number of energy zones for neutrinos of type n
!  unu(k)        : energy of energy zone k (MeV)
!  dunu(k)       : energy width of energy zone k (MeV)
!  psi0(j,k,n,ij_ray,ik_ray)
!                : zeroth moment of of the neutrino occupation probability
!                   for neutrinos of type n, energy zone k, radial zone j
!  psi1(j,k,n,ij_ray,ik_ray)
!                : first moment of of the neutrino occupation probability
!                   for neutrinos of type n, energy zone k, radial zone j
!  r_sphere(k,n,ij_ray,ik_ray)
!                : radius of the neutrinosphere for neutrinos of type n
!                   and energy k (cm).
!
!    Output arguments (common):
!
!   paf(1) = a0w     pafi(1) = da0w/di     pafpi(1,k) = da0w/dpsii(k)
!   paf(2) = b0w     pafi(2) = db0w/di     pafpi(2,k) = db0w/dpsii(k)
!   paf(3) = c0w     pafi(3) = dc0w/di     pafpi(3,k) = dc0w/dpsii(k)
!   paf(4) = a1w     pafi(4) = da1w/di     pafpi(4,k) = da1w/dpsii(k)
!   paf(5) = b1w     pafi(5) = db1w/di     pafpi(5,k) = db1w/dpsii(k)
!   paf(6) = c1w     pafi(6) = dc1w/di     pafpi(6,k) = dc1w/dpsii(k)
!
!    Subprograms called:
!      pairkrnl : interpolates the pair annihilation functions
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  edit_module, nu_dist_module, nu_energy_grid_module, pair_module,
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY : nez, nnu
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : kmev

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : unu, dunu, psi0, psi1, r_sphere
USE nu_energy_grid_module, ONLY : nnugp
USE pair_module, ONLY : paf, pafd, paft, pafy, pafp0, pafp1
USE prb_cntl_module, ONLY : ipair, rhopairemn, rhopairemx, rhopairtmn, rhopairtmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)               :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)               :: nx            ! radial array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rho    ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: t      ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: ye     ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: r      ! radius (cm)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                :: var_name

LOGICAL                           :: first_alloc = .true.
LOGICAL                           :: pair_off
LOGICAL                           :: first = .true.

INTEGER                           :: i             ! pair annihilation function index
INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! neutrino energy zone index

INTEGER                           :: kp            ! antineutrino energy zone index

INTEGER, DIMENSION(4)             :: nanti
INTEGER                           :: na            ! antineutrino flavor index

INTEGER                           :: istat         ! allocation status

REAL(KIND=double)                 :: enu           ! neutrino energy/kt
REAL(KIND=double)                 :: enubar        ! antineutrino energy/kt
REAL(KIND=double)                 :: tmev          ! temperature (MeV)

REAL(KIND=double)                 :: x1            ! psi1
REAL(KIND=double)                 :: x0a           ! psi1
REAL(KIND=double)                 :: x1a           ! psi1

REAL(KIND=double)                 :: c1            ! integrated scattering function
REAL(KIND=double)                 :: c1d           ! d(c1)/d(rho)
REAL(KIND=double)                 :: c1t           ! d(c1)/d(t)
REAL(KIND=double)                 :: c1y           ! d(c1)/d(ye)
REAL(KIND=double)                 :: gthp          ! integrated scattering function
REAL(KIND=double)                 :: gthpd         ! d(gthp)/d(rho)
REAL(KIND=double)                 :: gthpt         ! d(gthp)/d(t)
REAL(KIND=double)                 :: gthpy         ! d(gthp)/d(ye)
REAL(KIND=double)                 :: zlthp         ! integrated scattering function
REAL(KIND=double)                 :: zlthpd        ! d(zlthp)/d(rho)
REAL(KIND=double)                 :: zlthpt        ! d(zlthp)/d(t)
REAL(KIND=double)                 :: zlthpy        ! d(zlthp)/d(ye)
REAL(KIND=double)                 :: asthp         ! integrated scattering function
REAL(KIND=double)                 :: asthpd        ! d(asthp)/d(rho)
REAL(KIND=double)                 :: asthpt        ! d(asthp)/d(t)
REAL(KIND=double)                 :: asthpy        ! d(asthp)/d(ye)

REAL(KIND=double)                 :: arg           ! argument of an exponential
REAL(KIND=double)                 :: coefpp        ! coefficient for computing prodction from annihilation rates

REAL(KIND=double)                 :: f0p           ! zero moment of the pair-production function
REAL(KIND=double)                 :: f0pd          ! d(f0p)/d(rho)
REAL(KIND=double)                 :: f0pt          ! d(f0p)/d(t)
REAL(KIND=double)                 :: f0py          ! d(f0p)/d(ye)
REAL(KIND=double)                 :: f1p           ! first moment of the pair-production function
REAL(KIND=double)                 :: f1pd          ! d(f1p)/d(rho)
REAL(KIND=double)                 :: f1pt          ! d(f1p)/d(t)
REAL(KIND=double)                 :: f1py          ! d(f1p)/d(ye)

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

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: zlthpp0      ! d(zlthp)/d(psi0)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: c1p0         ! d(c1)/d(psi0)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: asthpp1      ! d(asthp)/d(psi1)
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: gthpp1       ! d(gthp)/d(psi0)

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

 1001 FORMAT (' Allocation problem for array ',a10,' in pairrate')
 2001 FORMAT (' Deallocation problem for array ',a10,' in pairrate')

!-----------------------------------------------------------------------        
!-----------------------------------------------------------------------        

IF ( istat /= 0 ) WRITE (nlog,1001)

!-----------------------------------------------------------------------        
!  Set pair annihilation functions to zero if
!     nnugp(n) =  0
!   or
!     ipair = 0
!
!  paf, pafd, paft, pafy, pafp0, pafp1 have been initialized to zero,
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0  .or.  ipair == 0 ) RETURN

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (w2dw(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w2dw      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (zlthpp0(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zlthpp0   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (c1p0(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c1p0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (asthpp1(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'asthpp1   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (gthpp1(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gthpp1    '; WRITE (nlog,1001) var_name; END IF

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

zlthpp0            = zero
c1p0               = zero
asthpp1            = zero
gthpp1             = zero

f0a                = zero
f0ad               = zero
f0at               = zero
f0ay               = zero
f1a                = zero
f1ad               = zero
f1at               = zero
f1ay               = zero

!-----------------------------------------------------------------------
!                      ||||| Loop over j |||||
!-----------------------------------------------------------------------

outer: DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Set pair annihilation functions to zero if rho outside specified
!   boundaries.
!-----------------------------------------------------------------------

  pair_off         = .false.
 
  IF ( n <  3  .and.  rho(j) < rhopairemn ) pair_off = .true.
  IF ( n <  3  .and.  rho(j) > rhopairemx ) pair_off = .true.
  IF ( n >= 3  .and.  rho(j) < rhopairtmn ) pair_off = .true.
  IF ( n >= 3  .and.  rho(j) > rhopairtmx ) pair_off = .true.

  IF ( pair_off ) THEN
    DO k = 1,nez
      DO i = 1,6
        paf(i,j,k,n)        = zero
        pafd(i,j,k,n)       = zero
        paft(i,j,k,n)       = zero
        pafy(i,j,k,n)       = zero
        DO kp = 1,nez
          pafp0(i,kp,j,k,n) = zero
          pafp1(i,kp,j,k,n) = zero
        END DO ! kp
      END DO ! i
    END DO ! k
    CYCLE ! j loop
  END IF ! pair_off

!-----------------------------------------------------------------------
!                      ||||| Loop over k |||||
!-----------------------------------------------------------------------

  DO k = 1,nnugp(n)

!-----------------------------------------------------------------------
!  Compute variables related to n-type neutrinos.
!
!  tmev       : t (MeV)
!  enu        : neutrino energy (MeV)
!  enubar     : antineutrino energy (MeV)
!  na         : neutrino index corresponding to the antiparticle of
!                the neutrino corresponding to index n
!-----------------------------------------------------------------------

    enu            = unu(j,k)
    tmev           = t(j) * kmev
    na             = nanti(n)

!-----------------------------------------------------------------------
!  Initialize electron-positron pair annihilation functions
!-----------------------------------------------------------------------

    zlthp          = zero
    asthp          = zero
    c1             = zero
    gthp           = zero
    zlthpd         = zero
    asthpd         = zero
    c1d            = zero
    gthpd          = zero
    zlthpt         = zero
    asthpt         = zero
    c1t            = zero
    gthpt          = zero
    zlthpy         = zero
    asthpy         = zero
    c1y            = zero
    gthpy          = zero

!-----------------------------------------------------------------------
!  Get electron-positron pair annihilation kernels
!-----------------------------------------------------------------------

    CALL pairkrnl(n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j), f0a, f1a, &
& f0ad, f1ad, f0at, f1at, f0ay, f1ay, nez )

    x1             = half * ( psi1(j-1,k,n ) + psi1(j,k,n ) )

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kp = 1,nnugp(n)

      w2dw(kp)     = unu(j,kp)**2 * dunu(j,kp)
      enubar       = unu(j,kp)
      x0a          = psi0(j,kp,na)
      x1a          = half * ( psi1(j-1,kp,na) + psi1(j,kp,na) )

      arg          = ( enu + enubar )/tmev
      coefpp       = fexp( -arg )

!-----------------------------------------------------------------------
!  Compute electron-positron pair production kernals.
!
!  f0p and f1p have dimensions of 1 /[ energy**3 length ].
!  zt, asp, c1, gp, and the scf(i) have dimensions of 1 /[ length ].
!----------------------------------------------------------------------

      f0p          = coefpp * f0a(kp)
      f1p          = coefpp * f1a(kp)
      f0pd         = coefpp * f0ad(kp)
      f1pd         = coefpp * f1ad(kp)
      f0pt         = coefpp * ( f0at(kp) + f0a(kp) * ( arg/t(j) ) )
      f1pt         = coefpp * ( f1at(kp) + f1a(kp) * ( arg/t(j) ) )
      f0py         = coefpp * f0ay(kp)
      f1py         = coefpp * f1ay(kp)

!-----------------------------------------------------------------------
!  Compute angular correction factors to electron-positron pair
!   annihilation kernals.
!-----------------------------------------------------------------------

      psib10       = DMIN1( DABS( x1a)/( x0a + epsilon ), one )
      psi10        = DMIN1( DABS( x1 )/( psi0(j,k,n ) + epsilon ), one )
      ac1          = one - psi10 * psib10/6.d0
      rsphdr       = DMIN1( r_sphere(k ,n ,ij_ray,ik_ray)/r(j), rsinmin )
      z            = DSQRT( one - rsphdr**2  )
      rsphbdr      = DMIN1( r_sphere(kp,na,ij_ray,ik_ray)/r(j), rsinmin )
      zb           = DSQRT( one - rsphbdr**2 )
      zav          = ( one + z  )/2.d0
      zbav         = ( one + zb )/2.d0
      z2av         = ( one + z  + z**2  )/3.d0
      zb2av        = ( one + zb + zb**2 )/3.d0
      ac2          = 0.75d0 * ( one - 2.d0 * zav * zbav + z2av * zb2av + 0.5d0  * ( one - z2av ) * ( one - zb2av ) )
      ac           = DMAX1( ac1, ac2, zero )

!-----------------------------------------------------------------------
!  Correct electron-positron pair annihilation kernals
!-----------------------------------------------------------------------

      f0a(kp)      = f0a(kp)  * ac
      f1a(kp)      = f1a(kp)  * ac
      f0ad(kp)     = f0ad(kp) * ac
      f1ad(kp)     = f1ad(kp) * ac
      f0at(kp)     = f0at(kp) * ac
      f1at(kp)     = f1at(kp) * ac
      f0ay(kp)     = f0ay(kp) * ac
      f1ay(kp)     = f1ay(kp) * ac

!-----------------------------------------------------------------------
!  Compute electron-positron pair production functions.
!
!  zlthp, asthp, c1, and gthp have dimensions of 1 /[ length ].
!-----------------------------------------------------------------------

      zlthp        = zlthp  + ( f0a(kp) * x0a + f0p * ( one - x0a ) )   * w2dw(kp)
      asthp        = asthp  + x1a * ( f1a(kp)  - f1p  )                 * w2dw(kp)
      c1           = c1     + f0p * ( one - x0a )                       * w2dw(kp)
      gthp         = gthp   - x1a * f1p                                 * w2dw(kp)

      zlthpd       = zlthpd + ( f0ad(kp) * x0a + f0pd * ( one - x0a ) ) * w2dw(kp)
      asthpd       = asthpd + x1a  * ( f1ad(kp) - f1pd )                * w2dw(kp)
      c1d          = c1d    + f0pd * ( one - x0a )                      * w2dw(kp)
      gthpd        = gthpd  - x1a  * f1pd                               * w2dw(kp)

      zlthpt       = zlthpt + ( f0at(kp) * x0a + f0pt * ( one - x0a ) ) * w2dw(kp)
      asthpt       = asthpt + x1a  * ( f1at(kp) - f1pt )                * w2dw(kp)
      c1t          = c1t    + f0pt * ( one - x0a )                      * w2dw(kp)
      gthpt        = gthpt  - x1a  * f1pt                               * w2dw(kp)

      zlthpy       = zlthpy + ( f0ay(kp) * x0a + f0py * ( one - x0a ) ) * w2dw(kp)
      asthpy       = asthpy + x1a  * ( f1ay(kp) - f1py )                * w2dw(kp)
      c1y          = c1y    + f0py * ( one - x0a )                      * w2dw(kp)
      gthpy        = gthpy  - x1a  * f1py                               * w2dw(kp)

      zlthpp0(kp)  = ( f0a(kp) - f0p )                                  * w2dw(kp)
      asthpp1(kp)  = ( f1a(kp) - f1p )                                  * w2dw(kp)
      c1p0(kp)     = - f0p                                              * w2dw(kp)
      gthpp1(kp)   = - f1p                                              * w2dw(kp)

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO ! kp

!-----------------------------------------------------------------------
!  Compute a0w, b0w, c0w, a1w, b1w, c1w, and derivatives
!
!   paf(1) = a0w     pafi(1) = da0w/di     pafpi(1,k) = da0w/dpsii(k)
!   paf(2) = b0w     pafi(2) = db0w/di     pafpi(2,k) = db0w/dpsii(k)
!   paf(3) = c0w     pafi(3) = dc0w/di     pafpi(3,k) = dc0w/dpsii(k)
!   paf(4) = a1w     pafi(4) = da1w/di     pafpi(4,k) = da1w/dpsii(k)
!   paf(5) = b1w     pafi(5) = db1w/di     pafpi(5,k) = db1w/dpsii(k)
!   paf(6) = c1w     pafi(6) = dc1w/di     pafpi(6,k) = dc1w/dpsii(k)
!
!  a0w, etc have dimensions of 1 /[ length ].
!-----------------------------------------------------------------------

    paf(1,j,k,n)         = -zlthp
    paf(2,j,k,n)         = -asthp/3.d0
    paf(3,j,k,n)         =  c1
    paf(4,j,k,n)         = -asthp
    paf(5,j,k,n)         = -zlthp
    paf(6,j,k,n)         =  gthp

    pafd(1,j,k,n)        = -zlthpd
    pafd(2,j,k,n)        = -asthpd/3.d0
    pafd(3,j,k,n)        =  c1d
    pafd(4,j,k,n)        = -asthpd
    pafd(5,j,k,n)        = -zlthpd
    pafd(6,j,k,n)        =  gthpd

    paft(1,j,k,n)        = -zlthpt
    paft(2,j,k,n)        = -asthpt/3.d0
    paft(3,j,k,n)        =  c1t
    paft(4,j,k,n)        = -asthpt
    paft(5,j,k,n)        = -zlthpt
    paft(6,j,k,n)        =  gthpt

    pafy(1,j,k,n)        = -zlthpy
    pafy(2,j,k,n)        = -asthpy/3.d0
    pafy(3,j,k,n)        =  c1y
    pafy(4,j,k,n)        = -asthpy
    pafy(5,j,k,n)        = -zlthpy
    pafy(6,j,k,n)        =  gthpy

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kp = 1,nnugp(n)

      pafp0(1,kp,j,k,n)  = -zlthpp0(kp)
      pafp0(2,kp,j,k,n)  =  zero
      pafp0(3,kp,j,k,n)  =  c1p0(kp)
      pafp0(4,kp,j,k,n)  =  zero
      pafp0(5,kp,j,k,n)  = -zlthpp0(kp)
      pafp0(6,kp,j,k,n)  =  zero

      pafp1(1,kp,j,k,n)  =  zero
      pafp1(2,kp,j,k,n)  = -asthpp1(kp)/3.
      pafp1(3,kp,j,k,n)  =  zero
      pafp1(4,kp,j,k,n)  = -asthpp1(kp)
      pafp1(5,kp,j,k,n)  =  zero
      pafp1(6,kp,j,k,n)  =  gthpp1(kp)

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO ! kp

!-----------------------------------------------------------------------
!                 ||||| End loop over j and k |||||
!-----------------------------------------------------------------------

  END DO ! k
END DO outer

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (w2dw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w2dw      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (zlthpp0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zlthpp0   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (c1p0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c1p0      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (asthpp1, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'asthpp1   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (gthpp1, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'gthpp1    '; WRITE (nlog,2001) var_name; END IF

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
END SUBROUTINE pairrate

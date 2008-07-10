SUBROUTINE bremrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, r, &
& nx )
!-----------------------------------------------------------------------
!
!    File:         bremrate
!    Module:       bremrate
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/21/97
!
!    Purpose:
!      The zero and first legendre coefs for the n-type neutrino bremsstrahlung
!       pair annihilation and production functions are computed here.
!      These are included in the multi-group diffusion equations, which have
!       the form
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
!      Neutrino bremsstrahlung pair annihilation and production contributes to the terms a0w,
!       a1w, b0w, b1w, c0w, c1w as follows:
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
!      The Legendre moments of the neutrino bremsstrahlung pair annihilation
!       and production functions are given by
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
!      The integrations are performed in subroutine bremcal.
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
!  ibrem         : 0, neutrino bremsstrahlung pair annihilation and production turned off;
!                   bremsstrahlung pair annihilation subroutines are bypassed; bremsstrahlung
!                   pair annihilation function arrays, if used, must be zeroed elsewhere
!                   1, neutrino bremsstrahlung pair annihilation and production included
!  rhobrememn    : rho < rhopairemn: e-neutrino-antineutrino bremsstrahlung pair annihilation
!                   and production omitted
!  rhobrememx    : rho > rhopairemx: e-neutrino-antineutrino bremsstrahlung pair annihilation
!                   and production omitted
!  rhobremtmn    : rho < rhopairtmn: x-neutrino-antineutrino bremsstrahlung pair annihilation
!                   and production omitted
!  rhobremtmx    : rho > rhopairtmx: x-neutrino-antineutrino bremsstrahlung pair annihilation
!                   and production omitted
!  nnugp(n)      : number of energy zones for neutrinos of type n
!  unu(k)        : energy of energy zone k (MeV)
!  dunu(k)       : energy width of energy zone k (MeV)
!  psi0(j,k,n)   : zeroth moment of of the neutrino occupation probability for neutrinos of type n,
!                   energy zone k, radial zone j
!  psi1(j,k,n)   : first moment of of the neutrino occupation probability for neutrinos of type n,
!                   energy zone k, radial zone j
!  r_sphere(k,n) : radius of the neutrinosphere for neutrinos of type n and energy k (cm).
!
!    Output arguments (common):
!
!   baf(1) = a0w     bafi(1) = da0w/di     bafpi(1,k) = da0w/dpsii(k)
!   baf(2) = b0w     bafi(2) = db0w/di     bafpi(2,k) = db0w/dpsii(k)
!   baf(3) = c0w     bafi(3) = dc0w/di     bafpi(3,k) = dc0w/dpsii(k)
!   baf(4) = a1w     bafi(4) = da1w/di     bafpi(4,k) = da1w/dpsii(k)
!   baf(5) = b1w     bafi(5) = db1w/di     bafpi(5,k) = db1w/dpsii(k)
!   baf(6) = c1w     bafi(6) = dc1w/di     bafpi(6,k) = dc1w/dpsii(k)
!
!    Subprograms called:
!      bremkrnl
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  brem_module, edit_module, nu_dist_module, nu_energy_grid_module,
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nez
USE numerical_module, ONLY : zero, half, one, epsilon
USE physcnst_module, ONLY : kmev

USE brem_module, ONLY : baf, bafd, baft, bafy, bafp0, bafp1
USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : unu, dunu, psi0, psi1, r_sphere
USE nu_energy_grid_module, ONLY : nnugp
USE prb_cntl_module, ONLY : ibrem, rhobrememn, rhobrememx, rhobremtmn, rhobremtmx

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

LOGICAL                           :: brem_off
LOGICAL                           :: first = .true.

INTEGER                           :: i             ! pair annihilation function index
INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! neutrino energy zone index

INTEGER                           :: kb            ! antineutrino energy zone index

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

 1001 FORMAT (' Allocation problem for array ',a10,' in bremrate')
 2001 FORMAT (' Deallocation problem for array ',a10,' in bremrate')

!-----------------------------------------------------------------------        
!-----------------------------------------------------------------------        

!-----------------------------------------------------------------------
!  Set bremsstrahlung pair annihilation functions to zero if
!     nnugp(n) =  0
!   or
!     ibrem = 0
!
!  baf, bafd, baft, bafy, bafp0, and bafp1 have been initialized to zero.
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0  .or.  ibrem == 0 ) RETURN

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

w2dw                        = zero

zlthpp0                     = zero
c1p0                        = zero
asthpp1                     = zero
gthpp1                      = zero

f0a                         = zero
f0ad                        = zero
f0at                        = zero
f0ay                        = zero
f1a                         = zero
f1ad                        = zero
f1at                        = zero
f1ay                        = zero

!-----------------------------------------------------------------------
!                      ||||| Loop over j |||||
!-----------------------------------------------------------------------

DO j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Set bremsstrahlung pair annihilation functions to zero if rho outside
!   specified boundaries.
!-----------------------------------------------------------------------

  brem_off                  = .false.

  IF ( n <  3  .and.  rho(j) < rhobrememn ) brem_off = .true.
  IF ( n <  3  .and.  rho(j) > rhobrememx ) brem_off = .true.
  IF ( n >= 3  .and.  rho(j) < rhobremtmn ) brem_off = .true.
  IF ( n >= 3  .and.  rho(j) > rhobremtmx ) brem_off = .true.

  IF ( brem_off ) THEN
    DO k = 1,nez
      DO i = 1,6
        baf(i,j,k,n)        = zero
        bafd(i,j,k,n)       = zero
        baft(i,j,k,n)       = zero
        bafy(i,j,k,n)       = zero
        DO kb = 1,nez
          bafp0(i,kb,j,k,n) = zero
          bafp1(i,kb,j,k,n) = zero
        END DO
      END DO
    END DO
    CYCLE
  END IF ! brem_off

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
!  Initialize bremsstrahlung pair annihilation functions
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
!  Get bremsstrahlung pair annihilation kernels
!-----------------------------------------------------------------------

    CALL bremkrnl( n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j), f0a, &
&    f1a, f0ad, f1ad, f0at, f1at, f0ay, f1ay, nez )

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kb = 1,nnugp(n)

      w2dw(kb)     = unu(j,kb)**2 * dunu(j,kb)
      enubar       = unu(j,kb)
      x0a          = psi0(j,kb,na)
      x1a          = half * ( psi1(j-1,kb,na) + psi1(j,kb,na) )
      x1           = half * ( psi1(j-1,kb,n ) + psi1(j,kb,n ) )

      arg          = ( enu + enubar )/tmev
      coefpp       = fexp( -arg )

!-----------------------------------------------------------------------
!  Compute bremsstrahlung pair production kernals.
!
!  f0p and f1p have dimensions of 1 /[ energy**3 length ].
!  zt, asp, c1, gp, and the scf(i) have dimensions of 1 /[ length ].
!-----------------------------------------------------------------------

      f0p          = coefpp * f0a(kb)
      f1p          = coefpp * f1a(kb)
      f0pd         = coefpp * f0ad(kb)
      f1pd         = coefpp * f1ad(kb)
      f0pt         = coefpp * ( f0at(kb) + f0a(kb) * ( arg/t(j) ) )
      f1pt         = coefpp * ( f1at(kb) + f1a(kb) * ( arg/t(j) ) )
      f0py         = coefpp * f0ay(kb)
      f1py         = coefpp * f1ay(kb)

!-----------------------------------------------------------------------
!  Compute angular correction factors to bremsstrahlung pair
!   annihilation kernals.
!-----------------------------------------------------------------------


      psib10       = DMIN1( dabs( x1a)/( x0a + epsilon ), one )
      psi10        = DMIN1( dabs( x1 )/( psi0(j,k,n ) + epsilon ), one )
      ac1          = one - psi10 * psib10/6.d0
      rsphdr       = DMIN1( r_sphere(k ,n ,ij_ray,ik_ray)/r(j), rsinmin )
      z            = DSQRT( one - rsphdr**2  )
      rsphbdr      = DMIN1( r_sphere(kb,na,ij_ray,ik_ray)/r(j), rsinmin )
      zb           = dsqrt( one - rsphbdr**2 )
      zav          = ( one + z  )/2.d0
      zbav         = ( one + zb )/2.d0
      z2av         = ( one + z  + z**2  )/3.d0
      zb2av        = ( one + zb + zb**2 )/3.d0
      ac2          = 0.75d0 * ( one - 2.d0 * zav * zbav + z2av * zb2av + 0.5d0  * ( one - z2av ) * ( one - zb2av ) )
      ac           = DMAX1( ac1, ac2, zero )

!-----------------------------------------------------------------------
!  Correct bremsstrahlung pair annihilation kernals
!-----------------------------------------------------------------------

      f0a(kb)      = f0a(kb)  * ac
      f1a(kb)      = f1a(kb)  * ac
      f0ad(kb)     = f0ad(kb) * ac
      f1ad(kb)     = f1ad(kb) * ac
      f0at(kb)     = f0at(kb) * ac
      f1at(kb)     = f1at(kb) * ac
      f0ay(kb)     = f0ay(kb) * ac
      f1ay(kb)     = f1ay(kb) * ac

!-----------------------------------------------------------------------
!  Compute bremsstrahlung pair production functions.
!
!  zlthp, asthp, c1, and gthp have dimensions of 1 /[ length ].
!-----------------------------------------------------------------------

      zlthp        = zlthp  + ( f0a(kb) * x0a + f0p * ( one - x0a ) )   * w2dw(kb)
      asthp        = asthp  + x1a * ( f1a(kb)  - f1p  )                 * w2dw(kb)
      c1           = c1     + f0p * ( one - x0a )                       * w2dw(kb)
      gthp         = gthp   - x1a * f1p                                 * w2dw(kb)

      zlthpd       = zlthpd + ( f0ad(kb) * x0a + f0pd * ( one - x0a ) ) * w2dw(kb)
      asthpd       = asthpd + x1a  * ( f1ad(kb) - f1pd )                * w2dw(kb)
      c1d          = c1d    + f0pd * ( one - x0a )                      * w2dw(kb)
      gthpd        = gthpd  - x1a  * f1pd                               * w2dw(kb)

      zlthpt       = zlthpt + ( f0at(kb) * x0a + f0pt * ( one - x0a ) ) * w2dw(kb)
      asthpt       = asthpt + x1a  * ( f1at(kb) - f1pt )                * w2dw(kb)
      c1t          = c1t    + f0pt * ( one - x0a )                      * w2dw(kb)
      gthpt        = gthpt  - x1a  * f1pt                               * w2dw(kb)

      zlthpy       = zlthpy + ( f0ay(kb) * x0a + f0py * ( one - x0a ) ) * w2dw(kb)
      asthpy       = asthpy + x1a  * ( f1ay(kb) - f1py )                * w2dw(kb)
      c1y          = c1y    + f0py * ( one - x0a )                      * w2dw(kb)
      gthpy        = gthpy  - x1a  * f1py                               * w2dw(kb)

      zlthpp0(kb)  = ( f0a(kb) - f0p )                                  * w2dw(kb)
      asthpp1(kb)  = ( f1a(kb) - f1p )                                  * w2dw(kb)
      c1p0(kb)     = - f0p                                              * w2dw(kb)
      gthpp1(kb)   = - f1p                                              * w2dw(kb)

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO

!-----------------------------------------------------------------------
!  Compute a0w, b0w, c0w, a1w, b1w, c1w, and derivatives
!
!   baf(1) = a0w     pafi(1) = da0w/di     pafpi(1,k) = da0w/dpsii(k)
!   baf(2) = b0w     pafi(2) = db0w/di     pafpi(2,k) = db0w/dpsii(k)
!   baf(3) = c0w     pafi(3) = dc0w/di     pafpi(3,k) = dc0w/dpsii(k)
!   baf(4) = a1w     pafi(4) = da1w/di     pafpi(4,k) = da1w/dpsii(k)
!   baf(5) = b1w     pafi(5) = db1w/di     pafpi(5,k) = db1w/dpsii(k)
!   baf(6) = c1w     pafi(6) = dc1w/di     pafpi(6,k) = dc1w/dpsii(k)
!
!  a0w, etc have dimensions of 1 /[ length ].
!-----------------------------------------------------------------------

    baf(1,j,k,n)             = -zlthp
    baf(2,j,k,n)             = -asthp/3.d0
    baf(3,j,k,n)             =  c1
    baf(4,j,k,n)             = -asthp
    baf(5,j,k,n)             = -zlthp
    baf(6,j,k,n)             =  gthp

    bafd(1,j,k,n)            = -zlthpd
    bafd(2,j,k,n)            = -asthpd/3.d0
    bafd(3,j,k,n)            =  c1d
    bafd(4,j,k,n)            = -asthpd
    bafd(5,j,k,n)            = -zlthpd
    bafd(6,j,k,n)            =  gthpd

    baft(1,j,k,n)            = -zlthpt
    baft(2,j,k,n)            = -asthpt/3.d0
    baft(3,j,k,n)            =  c1t
    baft(4,j,k,n)            = -asthpt
    baft(5,j,k,n)            = -zlthpt
    baft(6,j,k,n)            =  gthpt

    bafy(1,j,k,n)            = -zlthpy
    bafy(2,j,k,n)            = -asthpy/3.d0
    bafy(3,j,k,n)            =  c1y
    bafy(4,j,k,n)            = -asthpy
    bafy(5,j,k,n)            = -zlthpy
    bafy(6,j,k,n)            =  gthpy

!-----------------------------------------------------------------------
!                      ||||| Loop over kp |||||
!-----------------------------------------------------------------------

    DO kb = 1,nnugp(n)

      bafp0(1,kb,j,k,n)      = -zlthpp0(kb)
      bafp0(2,kb,j,k,n)      =  zero
      bafp0(3,kb,j,k,n)      =  c1p0(kb)
      bafp0(4,kb,j,k,n)      =  zero
      bafp0(5,kb,j,k,n)      = -zlthpp0(kb)
      bafp0(6,kb,j,k,n)      =  zero

      bafp1(1,kb,j,k,n)      =  zero
      bafp1(2,kb,j,k,n)      = -asthpp1(kb)/3.
      bafp1(3,kb,j,k,n)      =  zero
      bafp1(4,kb,j,k,n)      = -asthpp1(kb)
      bafp1(5,kb,j,k,n)      =  zero
      bafp1(6,kb,j,k,n)      =  gthpp1(kb)

!-----------------------------------------------------------------------
!                    ||||| End loop over kp |||||
!-----------------------------------------------------------------------

    END DO

!-----------------------------------------------------------------------
!                 ||||| End loop over j and k |||||
!-----------------------------------------------------------------------

  END DO
END DO

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
  IF ( istat /= 0 ) THEN; var_name = 'f1a      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1ad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ad     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1at, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1at     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f1ay, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f1ay     '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE bremrate

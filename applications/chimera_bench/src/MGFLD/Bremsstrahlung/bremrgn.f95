SUBROUTINE bremrgn( j, ij_ray, ik_ray, idd, itt, iyy, ida, ita, iya )
!-----------------------------------------------------------------------
!
!    File:         bremrgn
!    Module:       bremrgn
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/6/02
!
!    Purpose:
!      To set up the data for the computation of the neutrino-antineutrino
!       bremsstrahlung pair annihilation functions at the corners of the
!       "unit cell" in rho-T-ye state space, to call subroutine sctlgndv
!       for the computation of the scattering functions, and to load the
!       computed scattering functions in arrays for interpolation.
!
!-----------------------------------------------------------------------
!
!        The zero legendre coefs for the n-type neutrino-antineutrino
!         bremsstrahlung pair annihilation functions are computed here
!        These are included in the multi-group diffusion equations, which
!         have the form
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
!        Neutrino pair annihilation and production contributes to the terms
!         a0w, a1w, b0w, b1w, c0w, c1w as follows:
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
!        The integrations are performed in subroutine pairfnc, the quantities
!             [1 - Fe(e)]*[1 - F(w' + w - e)]*raLi(e,w,w')
!         and
!             [1 - Fe(e)]*[1 - F(w' + w - e)]*raLii(e,w,w')
!         are caculated in subroutines rjLi and rjLii, respectively.
!
!-----------------------------------------------------------------------
!
!    Subprograms called:
!        paircal
!
!    Input arguments:
!
!  j                     : radial zone number
!  ij_ray                : j-index of a radial ray
!  ik_ray                : k-index of a radial ray
!  idd                   : density grid index
!  itt                   : temperature grid index
!  iyy                   : electron fraction grid index
!  ida                   : bremsstrahlung pair annihilation function table
!                           density index
!  ita                   : bremsstrahlung pair annihilation function table
!                           temperature index
!  iya                   : bremsstrahlung pair annihilation function table
!                           electron fraction index
!
!    Input arguments (common):
!
!  rhobrememn            : density below which n_e - n_ebar bremsstrahlung
!                           pair annihilation is not computed (pair annihilation
!                           functions are zeroed).
!  rhobrememn            : density above which n_e - n_ebar bremsstrahlung pair
!                           annihilation is not computed (pair annihilation
!                           functions are zeroed).
!  rhobremtmn            : density below which n_x - n_xbar bremsstrahlung pair
!                           annihilation is not computed (pair annihilation
!                           functions are zeroed).
!  rhobremtmn            : density above which n_x - n_xbar bremsstrahlung pair
!                           annihilation is not computed (pair annihilation
!                           functions are zeroed).
!  nnugp(n)              : number of energy zones for neutrinos of type n
!  unu(k)                : energy of energy zone k (MeV)
!  dunu(k)               : energy width of energy zone k (MeV)
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and
!                          ye = 0 for zone j
!  idty(j,ij_ray,ik_ray) : index for dgrid, tgrid, and ygrid for zone j
!  estble(i,j,ida,ita,iya,ij_ray,ik_ray)
!                        : equation of state table array
!  brema0                : pair annihilation kernals
!
!    Include files:
!  kind_module, array_module, numerical_module
!  brem_module, eos_snc_x_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : ij_ray_dim
USE numerical_module, ONLY : one, epsilon

USE brem_module, ONLY : brema0
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, estble, escnst, idty
USE nu_dist_module, ONLY : unu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

INTEGER, INTENT(in)              :: idd           ! density do index
INTEGER, INTENT(in)              :: itt           ! temperature do index
INTEGER, INTENT(in)              :: iyy           ! electron fraction do index

INTEGER, INTENT(in)              :: ida           ! density cube index
INTEGER, INTENT(in)              :: ita           ! temperature cube index
INTEGER, INTENT(in)              :: iya           ! electron fraction cube index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: k             ! neutrino energy zone index
INTEGER                          :: kb            ! antineutrino energy zone index

INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

REAL(KIND=double)                :: rho           ! density (g/cm**3)
REAL(KIND=double)                :: t             ! temperature (K)
REAL(KIND=double)                :: ye            ! electron fraction

REAL(KIND=double)                :: xn            ! free neutron mass fraction
REAL(KIND=double)                :: xp            ! free proton mass fraction
REAL(KIND=double)                :: enu           ! incoming neutrino energy
REAL(KIND=double)                :: enubar        ! incoming antineutrino energy
REAL(KIND=double)                :: s_a           ! bremsstrahlung pair annihilation kernel

REAL(KIND=double)                :: f10           ! function 10**

EXTERNAL f10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Combine indices to uniquely index the radial ray
!-----------------------------------------------------------------------

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!  Compute independent variables at grid point
!-----------------------------------------------------------------------

rho                = f10   ( DBLE(idd)/dgrid(idty(j,ij_ray,ik_ray)) )
t                  = f10   ( DBLE(itt)/tgrid(idty(j,ij_ray,ik_ray)) )
ye                 = one - ( DBLE(iyy)/ygrid(idty(j,ij_ray,ik_ray)) )

!-----------------------------------------------------------------------
!  Fetch equation of state quantities needed for the computation of the
!   bremsstrahlung pair annihilation functions.
!-----------------------------------------------------------------------

xn                 = f10( DBLE( estble(7 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 7,j,ij_ray,ik_ray) - epsilon
xp                 = f10( DBLE( estble(8 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 8,j,ij_ray,ik_ray) - epsilon

!-----------------------------------------------------------------------
!  Neutrino-antineutrino bremsstrahlung pair annihilation functions.
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  enu              = unu(j,k)
  DO kb = 1,nnugpmx
    enubar         = unu(j,kb)
    CALL bremcal( enu, enubar, rho, t, xn, xp, s_a )
    brema0(j,k,kb,ida,ita,iya,i_ray)  = real(dlog10( dabs(s_a ) + epsilon ))
  END DO
END DO

RETURN
END SUBROUTINE bremrgn

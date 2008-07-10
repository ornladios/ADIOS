SUBROUTINE scatnnrgn( j, ij_ray, ik_ray, idd, itt, iyy, ida, ita, iya )
!-----------------------------------------------------------------------
!
!    File:         scatnnrgn
!    Module:       scatnnrgn
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/10/03
!
!    Purpose:
!      To set up the data for the computation of the neutrino- nucleon
!       inelastic scattering functions at the corners of the "unit cell"
!       in rho-T-ye state space, to call subroutine sctnncal for the
!       computation of the scattering functions, and to load the computed
!       scattering functions in arrays for interpolation.
!
!-----------------------------------------------------------------------
!
!      The zero and first legendre coefs for the n-type neutrino scattering
!       functions are computed here
!      These are included in the multi-group diffusion equations, which
!       have the form
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
!      Neutrino scattering contributes to the terms a0w,a1w, b0w, b1w,
!       c0w, c1w as follows:
!
!            a0w  = -K   Int{w2'dw'[ sctin0(w,w')psi0(w') + sctot0(w,w')( 1 - psi0(w') ) ] }
!            b0w  = -K/3 Int{w2'dw'[ sctin1(w,w') - sctot1(w,w') ]psi1(w') }
!            c0w  =  K   Int{w2'dw'[ sctin0(w,w')psi0(w') ] }
!            a1w  = -K   Int{w2'dw'[ sctin1(w,w') - sctot1(w,w') ]psi1(w') }
!            b1w  = -K   Int{w2'dw'[ sctin0(w,w')psi0(w') + sctot0(w,w')( 1 - psi0(w') ) ] }
!            c1w  =  K   Int{w2'dw'[ sctin1(w,w')psi1(w') ] }
!
!       where
!
!            K    = 2*pi/! * 1/(hc)**3
!
!       for neutrino-electron scattering.
!
!      The Legendre moments of the neutrino-electron scattering functions
!       are given by
!
!          phiLin = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 )
!            * cnes1(n)*hinLi(w,w') + cnes2(n)*hinLii(w,w')
!
!       where
!
!     hinLi(w,w') = Int{ de*Fe(e)*[1 - F(e + w - w')] * exp[(w - w')/kT]*hLi(e,w,w') }
!
!                 = (kT)**6 * Int{ d(e/kT) * Fe(e/kT)*[1 - F(e/kT + w/kT - w'/kT)]
!                  * exp[(w - w')/kT]*hLi(e/kT,w/kT,w'/kT) }
!
!      and hinLii(w,w') is defined likewise.
!
!         phiLout = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 ) * cnes1(n)*houtLi(w,w') + cnes2(n)*houtLii(w,w')
!
!      where
!
!    houtLi(w,w') = Int{ de*Fe(e)*[1 - F(e + w - w')] * hLo(e,w,w') }
!
!                 = (kT)**6 * Int{ d(e/kT) * Fe(e/kT)*[1 - F(e/kT + w/kT - w'/kT)]
!                  * hLo(e/kT,w/kT,w'/kT) }
!
!      and houtLii(w,w') is defined likewise.
!
!     The integrations are performed in subroutine sctnncal.
!
!-----------------------------------------------------------------------
!
!    Subprograms called:
!        sctnncal
!
!    Variables that must be passed through common:
!  nnugp(n)     : number of energy zones for neutrinos of type n
!  unu(k)       : energy of energy zone k (MeV)
!  dunu(k)      : energy width of energy zone k (MeV)
!  dgrid(idty(j,ij_ray,ik_ray))
!               : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray))
!               : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray))
!               : number of table entries in ye between ye = 0.5 and
!                  ye = 0 for zone j
! idty(j,ij_ray,ik_ray) : index for dgrid, tgrid, and ygrid for zone j
!  estble(i,j,ida,ita,iya,ij_ray,ik_ray)
!               : equation of state table array
!   sctnn0      : neutrino-nucleon inelastic scattering kernals
!
!    Input arguments:
!
!  j            : radial zone number
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  idd          : density grid index.
!  itt          : temperature grid index.
!  iyy          : electron fraction grid index.
!  ida          : density index of scattering function array
!                  sctnn0, i.e., indices of unit EOS table
!                  cube centered about state point.
!  ita          : temperature index of scattering function array
!                  sctnn0, i.e., indices of unit EOS table
!                  cube centered about state point.
!  iya          : electron fraction index of scattering function array
!                  sctnn0, i.e., indices of unit EOS table
!                  cube centered about state point.
!
!    Input arguments (common):
!
!  rhosctnnemn  : density below which e-neutrino-nucleon inelastic
!                  scattering is turned off.
!  rhosctnnemx  : density above which e-neutrino-nucleon inelastic
!                  scattering is turned off.
!  rhosctnntmn  : density below which t-neutrino-nucleon inelastic
!                  scattering is turned off.
!  rhosctnntmx  : density above which t-neutrino-nucleon inelastic
!                  scattering is turned off.
!
!    Include files:
!      kind_module, array_module, numerical_module
!      eos_snc_x_module, nu_dist_module, nu_energy_grid_module, scat_nn_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : ij_ray_dim
USE numerical_module, ONLY : one, epsilon

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, estble, escnst, idty
USE nu_dist_module, ONLY : unu, unub
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE scat_nn_module, ONLY : sctnn0

IMPLICIT NONE
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

INTEGER                          :: k             ! incomiong neutrino energy zone index
INTEGER                          :: kp            ! outgoing neutrino energy zone index

INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

REAL(KIND=double)                :: rho           ! density (g/cm**3)
REAL(KIND=double)                :: t             ! temperature (K)
REAL(KIND=double)                :: ye            ! electron fraction

REAL(KIND=double)                :: xn            ! neutron mass fraction
REAL(KIND=double)                :: xp            ! proton mass fraction

REAL(KIND=double)                :: enuin         ! incoming neutrino zone-centered energy
REAL(KIND=double)                :: enubk         ! incoming neutrino lower zone-edged energy
REAL(KIND=double)                :: enubkp1       ! incoming neutrino upper zone-edged energy
REAL(KIND=double)                :: enuout        ! outgoing neutrino zone-centered energy
REAL(KIND=double)                :: enubkp        ! outgoing neutrino lower zone-edged energy
REAL(KIND=double)                :: enubkpp1      ! outgoing neutrino upper zone-edged energy
REAL(KIND=double)                :: sct_nn        ! neutrino-neutron inelastic sct function

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
!   neutrino-nucleon inelastic scattering functions.
!-----------------------------------------------------------------------

xn                 = f10( DBLE( estble(7 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 7,j,ij_ray,ik_ray) - epsilon
xp                 = f10( DBLE( estble(8 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 8,j,ij_ray,ik_ray) - epsilon

!-----------------------------------------------------------------------
!  Neutrino-nucleon down and iso-e scattering kernals
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  enuin            = unu(j,k)
  enubk            = unub(j,k)
  enubkp1          = unub(j,k+1)

  DO kp = 1,k
    enuout         = unu(j,kp)
    enubkp         = unub(j,kp)
    enubkpp1       = unub(j,kp+1)

    CALL sctnncal( enuin, enuout, enubk, enubkp1, enubkp, enubkpp1, rho, &
&    t, xn, xp, sct_nn )
    sctnn0(j,k,kp,ida,ita,iya,i_ray) = REAL( LOG10( DABS(sct_nn) + epsilon ) )
    sctnn0(j,kp,k,ida,ita,iya,i_ray) = REAL( LOG10( DABS(sct_nn) + epsilon ) )

  END DO
END DO

RETURN
END SUBROUTINE scatnnrgn

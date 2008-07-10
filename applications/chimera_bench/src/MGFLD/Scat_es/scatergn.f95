SUBROUTINE scatergn( j, ij_ray, ik_ray, idd, itt, iyy, ida, ita, iya )
!-----------------------------------------------------------------------
!
!    File:         scatergn
!    Module:       scatergn
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/3/03
!
!    Purpose:
!      To set up the data for the computation of the neutrino-
!       electron scattering functions at the corners of the
!       "unit cell" in rho-T-ye state space, to call subroutine
!       sctlgndv for the computation of the scattering functions,
!       and to load the computed scattering functions in arrays
!       for interpolation.
!
!-----------------------------------------------------------------------
!
!      The zero and first legendre coefs for the n-type neutrino scattering
!       functions are computed here These are included in the multi-group
!       diffusion equations, which have the form
!
!            pis1 = -yt*( dpsi0/dr - a1w*psi0 - c1w)
!
!            dpsi0/dt/! + d(r2*psi1)/dr/3./r2
!                 = xw + yw*psi0 + zw*dpsi0/dr
!
!       where
!
!            1/yt = zltot = jw + yaw - b1w
!            xw   = jw + c0w + b0w*c1w*yt
!            yw   = -jw - yaw + a0w + b0w*a1w*yt
!            zw   = -b0w*yt
!
!      Neutrino scattering contributes to the terms a0w,a1w, b0w, b1w, c0w, c1w as follows:
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
!          phiLin = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 ) * cnes1(n)*hinLi(w,w') + cnes2(n)*hinLii(w,w')
!
!       where!
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
!    houtLi(w,w') = Int{ de*Fe(e)*[1 - F(e + w - w')]
!                  * hLo(e,w,w') }
!
!                 = (kT)**6 * Int{ d(e/kT) * Fe(e/kT)*[1 - F(e/kT + w/kT - w'/kT)] * hLo(e/kT,w/kT,w'/kT) }
!
!      and houtLii(w,w') is defined likewise.
!
!     The integrations are performed in subroutine sctgldnv.
!
!-----------------------------------------------------------------------
!
!    Subprograms called:
!        sctlgndv
!
!    Variables that must be passed through common:
!
!  nnugp(n)   : number of energy zones for neutrinos of type n
!  unu(k)     : energy of energy zone k (MeV)
!  dunu(k)    : energy width of energy zone k (MeV)
!  dgrid(idty(j,ij_ray,ik_ray))
!             : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray))
!             : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray))
!             : number of table entries in ye between ye = 0.5 and
!                ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)
!             : index for dgrid, tgrid, and ygrid for zone j
!  estble(i,j,ida,ita,iya,ij_ray,ik_ray)
!             : equation of state table array
!  scte0i,
!  scte0ii    : electron scattering kernals
!
!    Input arguments:
!
!  j          : radial zone number
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!  idd        : density grid index.
!  itt        : temperature grid index.
!  iyy        : electron fraction grid index.
!  ida        : density index of scattering function arrays scte0i and scte0ii,
!                i.e., indices of unit EOS table cube centered about state point.
!  ita        : temperature index of scattering function arrays scte0i and scte0ii,
!                i.e., indices of unit EOS table cube centered about state point.
!  iya        : electron fraction index of scattering function arrays scte0i and scte0ii,
!                i.e., indices of unit EOS table cube centered about state point.
!
!    Input arguments (common):
!
!  rhonesmn   : density below which neutrino electron scattering is not computed.
!  rhonesmx   : density above which neutrino electron scattering is not computed.
!
!    Include files:
!  kind_module,  array_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, nu_dist_module, nu_energy_grid_module,
!  scat_e_module
!
!-----------------------------------------------------------------------
!        Units
!
!  Gw         : G_fermi * ( mpc^{2} )^{2} dimensionless
!  Gw/mp**2   : MeV^{-2}
!  coc        : MeV^{-5} cm^{-1}
!  cxct       : MeV cm^{-1}
!  cxc        : MeV^{-1} cm^{-1}
!  cxc/wkp2   : MeV^{-3} cm^{-1}
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : ij_ray_dim
USE numerical_module, ONLY : one, epsilon
USE physcnst_module, ONLY : Gw, mp, pi, hbar, cvel, kmev

USE edit_module, ONLY : nprint
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, estble, escnst, idty
USE nu_dist_module, ONLY : unu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE scat_e_module, ONLY : scte0i, scte0ii

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

LOGICAL                          :: first = .true.

INTEGER                          :: k             ! incomiong neutrino energy zone index
INTEGER                          :: kp            ! outgoing neutrino energy zone index

INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

REAL(KIND=double)                :: rho           ! density (g/cm**3)
REAL(KIND=double)                :: t             ! temperature (K)
REAL(KIND=double)                :: ye            ! electron fraction

REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double)                :: cmpe          ! electron chemical potential
REAL(KIND=double)                :: eta           ! cmpe/tmev
REAL(KIND=double)                :: wk            ! incoming neutrino energy
REAL(KIND=double)                :: wk2           ! wk**2
REAL(KIND=double)                :: wkp2          ! outgoing neutrino energy squared
REAL(KIND=double)                :: enuin         ! incoming neutrino energy/tmev
REAL(KIND=double)                :: enuout        ! outgoing neutrino energy/tmev

REAL(KIND=double)                :: coc           ! 2*pi/!(2*pi*hbar*!)**3 !(hbar*!)**2*G**2/pi
REAL(KIND=double)                :: cxct          ! coc*t**6
REAL(KIND=double)                :: cxc           ! cxc/wk**2

REAL(KIND=double)                :: hin0i         ! zero moment of incoming scattering function type i
REAL(KIND=double)                :: hin0ii        ! zero moment of incoming scattering function type ii
REAL(KIND=double)                :: hout0i        ! zero moment of outgoing scattering function type i
REAL(KIND=double)                :: hout0ii       ! zero moment of outgoing scattering function type ii
REAL(KIND=double)                :: hin1i         ! first moment of incoming scattering function type i
REAL(KIND=double)                :: hin1ii        ! first moment of incoming scattering function type ii
REAL(KIND=double)                :: hout1i        ! first moment of outgoing scattering function type i
REAL(KIND=double)                :: hout1ii       ! first moment of outgoing scattering function type ii

REAL(KIND=double)                :: f10           ! function 10**

EXTERNAL f10

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  999 FORMAT (' nnugpmx = 0 in subroutine scatergn')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  coc              = 2.d0 * ( Gw/mp**2 )**2 * 1.d0/( 8.d0 * pi**3 * hbar * cvel )
  first            = .false.
END IF

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
!  Fetch equation of state quantities needed for the computation
!   of the neutrino scattering functions.
!-----------------------------------------------------------------------

cmpe              = f10( DBLE( estble(6 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                 - escnst( 6,j,ij_ray,ik_ray) - epsilon

!-----------------------------------------------------------------------
!  Compute equation of state quantities needed for the computation
!   of the neutrino scattering functions.
!
!  cxct has dimensions of [ energy / length ].
!-----------------------------------------------------------------------

tmev               = t * kmev
eta                = cmpe/tmev
cxct               = coc * (tmev)**6

!-----------------------------------------------------------------------
!  Neutrino-electron down and iso-e scattering kernals
!
!  cxc has dimensions of 1 /[ energy length ]
!  cxc/wje2 has dimensions of 1 /[ energy**3 length ]
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) THEN
  WRITE (nprint,999)
  STOP
END IF

DO k = 1,nnugpmx

  wk               = unu(j,k)
  wk2              = wk * wk
  cxc              = cxct/wk2
  enuin            = wk/tmev

  DO kp = 1,k

    enuout         = unu(j,kp)/tmev
      CALL sctlgndv( enuin, enuout, eta, hin0i, hin0ii, hout0i, hout0ii, &
&      hin1i, hin1ii, hout1i, hout1ii )

    wkp2           = unu(j,kp) * unu(j,kp)

    scte0i (j,k,kp,ida,ita,iya,i_ray) = REAL( log10( hout0i * cxc/wkp2 + epsilon ) )
    scte0ii(j,k,kp,ida,ita,iya,i_ray) = REAL( log10( hout0ii* cxc/wkp2 + epsilon ) )
    scte0i (j,kp,k,ida,ita,iya,i_ray) = REAL( log10( hout0i * cxc/wkp2 + epsilon ) )
    scte0ii(j,kp,k,ida,ita,iya,i_ray) = REAL( log10( hout0ii* cxc/wkp2 + epsilon ) )

  END DO
END DO

RETURN
END SUBROUTINE scatergn

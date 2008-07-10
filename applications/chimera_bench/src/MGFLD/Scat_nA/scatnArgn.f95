SUBROUTINE scatnArgn( j, ij_ray, ik_ray, idd, itt, iyy, ida, ita, iya )
!-----------------------------------------------------------------------
!
!    File:         scatnArgn
!    Module:       scatnArgn
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/10/03
!
!    Purpose:
!      To set up the data for the computation of the neutrino-
!       nucleus inelastic scattering functions at the corners of the
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
!      Neutrino scattering contributes to the terms a0w,a1w, b0w, b1w, c0w,
!       c1w as follows:
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
!    houtLi(w,w') = Int{ de*Fe(e)*[1 - F(e + w - w')] * hLo(e,w,w') }
!
!                 = (kT)**6 * Int{ d(e/kT) * Fe(e/kT)*[1 - F(e/kT + w/kT - w'/kT)]
!                  * hLo(e/kT,w/kT,w'/kT) }
!
!      and houtLii(w,w') is defined likewise.
!
!     The integrations are performed in subroutine sctnAcal.
!
!-----------------------------------------------------------------------
!
!    Subprograms called:
!        sctnAcal
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
!  sctnA0       : neutrino-nucleus inelastic scattering kernals
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
!                  sctnA0, i.e., indices of unit EOS table
!                  cube centered about state point.
!  ita          : temperature index of scattering function array
!                  sctnA0, i.e., indices of unit EOS table
!                  cube centered about state point.
!  iya          : electron fraction index of scattering function array
!                  sctnA0, i.e., indices of unit EOS table
!                  cube centered about state point.
!
!    Input arguments (common):
!
!  rhosctnAemn  : density below which e-neutrino-nucleus inelastic
!                  scattering is turned off.
!  rhosctnAemx  : density above which e-neutrino-nucleus inelastic
!                  scattering is turned off.
!  rhosctnAtmn  : density below which t-neutrino-nucleus inelastic
!                  scattering is turned off.
!  rhosctnAtmx  : density above which t-neutrino-nucleus inelastic
!                  scattering is turned off.
!
!    Include files:
!      kind_module, array_module, numerical_module, physcnst_module
!      edit_module, eos_snc_x_module, nu_dist_module, nu_energy_grid_module,
!      scat_a_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu, ij_ray_dim
USE numerical_module, ONLY : one, epsilon
USE physcnst_module, ONLY : Gw, mp, pi, hbar, cvel, kmev

USE edit_module, ONLY : nprint
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, estble, escnst, idty
USE nu_dist_module, ONLY : unu, unub
USE nu_energy_grid_module, ONLY : nnugp
USE scat_nA_module, ONLY : sctnA0

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
INTEGER                          :: n             ! neutrino flavor index

INTEGER                          :: nnugpmx       ! maximum neutrino group number

INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

REAL(KIND=double)                :: rho           ! density (g/cm**3)
REAL(KIND=double)                :: t             ! temperature (K)
REAL(KIND=double)                :: ye            ! electron fraction

REAL(KIND=double)                :: tmev          ! temperature (MeV)
REAL(KIND=double)                :: XH            ! heavy nuclei mass fraction
REAL(KIND=double)                :: A             ! mean nuclear mass number
REAL(KIND=double)                :: Z             ! mean nuclear charge number
REAL(KIND=double)                :: e_in          ! incoming neutrino energy/tmev
REAL(KIND=double)                :: e_in_l        ! incoming neutrino lower zone-edged energy
REAL(KIND=double)                :: e_in_u        ! incoming neutrino upper zone-edged energy
REAL(KIND=double)                :: e_out         ! outgoing neutrino energy/tmev
REAL(KIND=double)                :: e_out_l       ! outgoing neutrino lower zone-edged energy
REAL(KIND=double)                :: e_out_u       ! outgoing neutrino upper zone-edged energy
REAL(KIND=double)                :: phi_fm        ! neutrino-nucleus inelastic scattering function

REAL(KIND=double)                :: f10           ! function 10**

EXTERNAL f10

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  999 FORMAT (' nnugpmx = 0 in subroutine scatnArgn')

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
!   neutrino-nucleus inelastic scattering functions.
!-----------------------------------------------------------------------

XH                 = f10( DBLE( estble( 9,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst(  9,j,ij_ray,ik_ray) - epsilon
A                  = f10( DBLE( estble(10,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 10,j,ij_ray,ik_ray) - epsilon
Z                  = f10( DBLE( estble(11,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 11,j,ij_ray,ik_ray) - epsilon

!-----------------------------------------------------------------------
!  Compute equation of state quantities needed for the computation of
!   the neutrino-nucleus inelastic scattering functions.
!
!     cxct has dimensions of [ energy / length ].
!-----------------------------------------------------------------------

tmev               = t * kmev

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering kernals
!-----------------------------------------------------------------------

nnugpmx            = 0
DO n = 1,nnu
  nnugpmx          = MAX( nnugpmx, nnugp(n) )
END DO

IF ( nnugpmx == 0 ) THEN
  WRITE (nprint,999)
  STOP
END IF

!-----------------------------------------------------------------------
!
!        \\\\\ COMPUTE SCATTERING KERNELS AT CUBE CORNERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute dounscattering and isoenergetic scattering only for
!   out-scattering (k,kp) and in-scattering (kp,k).
!  Upscattering will be computed when needed from matter statistical
!   equilibrium relations.
!-----------------------------------------------------------------------


DO k = 1,nnugpmx

  e_in             = unu(j,k)
  e_in_l           = unub(j,k)
  e_in_u           = unub(j,k+1)

  DO kp = 1,nnugpmx

    e_out          = unu(j,kp)
    e_out_l        = unub(j,kp)
    e_out_u        = unub(j,kp+1)

    CALL sctnAcal( e_in, e_out, e_in_l, e_in_u, e_out_l, e_out_u, tmev, &
& XH, A, Z, rho, phi_fm)

    sctnA0 (j,k,kp,ida,ita,iya,i_ray) = REAL( log10( phi_fm + epsilon ) )

  END DO
END DO

RETURN
END SUBROUTINE scatnArgn

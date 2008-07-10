SUBROUTINE scatnrgn( j, ij_ray, ik_ray, idd, itt, iyy, ida, ita, iya )
!-----------------------------------------------------------------------
!
!    File:         scatnrgn
!    Module:       scatnrgn
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/09/03
!
!    Purpose:
!      To set up the data for the computation of the neutrino-
!       nucleon elastic scattering functions at the corners of the
!       "unit cell" in rho-T-ye state space, to call subroutine
!       sctncal for the computation of the scattering functions,
!       and to load the computed scattering functions in arrays
!       for interpolation.
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
!          phiLin = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 ) * cnes1(n)*hinLi(w,w') + cnes2(n)*hinLii(w,w')
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
!                 = (kT)**6 * Int{ d(e/kT)* Fe(e/kT)*[1 - F(e/kT + w/kT - w'/kT)]
!                  * hLo(e/kT,w/kT,w'/kT) }
!
!      and houtLii(w,w') is defined likewise.
!
!     The integrations are performed in subroutine sctgldnv.
!
!-----------------------------------------------------------------------
!
!    Subprograms called:
!        sctncal
!
!    Variables that must be passed through common:
!  nnugp(n)             : number of energy zones for neutrinos of type n
!  unu(k)               : energy of energy zone k (MeV)
!  dunu(k)              : energy width of energy zone k (MeV)
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and 
!                          ye = 0 for zone j
!  idty(j,ij_ray,ik_ray): index for dgrid, tgrid, and ygrid for zone j
!  estble(i,j,ida,ita,iya,ij_ray,ik_ray)
!                       : equation of state table array
!  sctn0                : neutrino-nucleon inelastic scattering kernals
!
!    Input arguments:
!
!  j                    : radial zone number
!  ij_ray               : j-index of a radial ray
!  ik_ray               : k-index of a radial ray
!  idd                  : density grid index.
!  itt                  : temperature grid index.
!  iyy                  : electron fraction grid index.
!  ida                  : density index of scattering function arrays sctn0, i.e.,
!                          indices of unit EOS table cube centered about state point.
!  ita                  : temperature index of scattering function array sctn0, i.e.,
!                          indices of unit EOS table cube centered about state point.
!  iya                  : electron fraction index of scattering function array sctn0, i.e.,
!                          indices of unit EOS table cube centered about state point.
!
!    Include files:
!      kind_module, array_module, numerical_module, physcnst_module
!      edit_module, eos_snc_x_module, nu_dist_module, nu_energy_grid_module,
!      scat_n_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : ij_ray_dim
USE numerical_module, ONLY : half, one, epsilon
USE physcnst_module, ONLY : ga, sin2W, kmev, dmnp, mn, mp

USE edit_module, ONLY : nprint
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, estble, escnst, idty
USE nu_dist_module, ONLY : unu, unub
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE scat_n_module, ONLY : sctn0, sctnb0

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
REAL(KIND=double)                :: cmpn          ! neutron chemical potential
REAL(KIND=double)                :: cmpp          ! proton chemical potential
REAL(KIND=double)                :: sm            ! target particle rest mass

REAL(KIND=double)                :: e_in          ! incoming neutrino zone-centered energy
REAL(KIND=double)                :: e_in_l        ! incoming neutrino lower zone-edged energy
REAL(KIND=double)                :: e_in_u        ! incoming neutrino upper zone-edged energy
REAL(KIND=double)                :: e_out         ! outgoing neutrino zone-centered energy
REAL(KIND=double)                :: e_out_l       ! outgoing neutrino lower zone-edged energy
REAL(KIND=double)                :: e_out_u       ! outgoing neutrino upper zone-edged energy

REAL(KIND=double)                :: cv_n          ! neutrino-neutron vector coupling constant
REAL(KIND=double)                :: ca_n          ! neutrino-neutron axial vector coupling constant
REAL(KIND=double)                :: cv_p          ! neutrino-proton vector coupling constant
REAL(KIND=double)                :: ca_p          ! neutrino-proton vector coupling constant

REAL(KIND=double)                :: sct0_n        ! zero momemt of the neutrino-neutron sct function
REAL(KIND=double)                :: sct1_n        ! first momemt of the neutrino-neutron sct function
REAL(KIND=double)                :: sct0_p        ! zero momemt of the neutrino-proton sct function
REAL(KIND=double)                :: sct1_p        ! first momemt of the neutrino-neutron sct function

REAL(KIND=double)                :: xi_p_wm       ! weak magnetism correction for neutrino-proton scattering
REAL(KIND=double)                :: xi_n_wm       ! weak magnetism correction for neutrino-neutron scattering
REAL(KIND=double)                :: xib_p_wm      ! weak magnetism correction for antineutrino-proton scattering
REAL(KIND=double)                :: xib_n_wm      ! weak magnetism correction for antineutrino-neutron scattering

REAL(KIND=double)                :: f10           ! function 10**

EXTERNAL f10

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  999 FORMAT (' nnugpmx = 0 in subroutine scatnrgn')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first) THEN
  cv_n             = - half
  ca_n             = - ga/2.d0
  cv_p             = half - 2.d0 * sin2W
  ca_p             = ga/2.d0
  first            = .false.
END IF

!-----------------------------------------------------------------------
!  Combine indices to uniquely index the radial ray
!-----------------------------------------------------------------------

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!
!             \\\\\ GET INPUT DATA AT CUBE CORNERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute independent variables at grid point
!-----------------------------------------------------------------------

rho                = f10   ( DBLE(idd)/dgrid(idty(j,ij_ray,ik_ray)) )
t                  = f10   ( DBLE(itt)/tgrid(idty(j,ij_ray,ik_ray)) )
ye                 = one - ( DBLE(iyy)/ygrid(idty(j,ij_ray,ik_ray)) )

tmev               = kmev * t

IF ( nnugpmx == 0 ) THEN
  WRITE (nprint,999)
  STOP
END IF

!-----------------------------------------------------------------------
!  Fetch equation of state quantities needed for the computation of the
!   neutrino-nucleon elastic scattering functions.
!-----------------------------------------------------------------------

cmpn               = f10( DBLE( estble(4 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 4,j,ij_ray,ik_ray) - epsilon + dmnp + mn
cmpp               = f10( DBLE( estble(5 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 5,j,ij_ray,ik_ray) - epsilon + dmnp + mp

!-----------------------------------------------------------------------
!
!        \\\\\ COMPUTE SCATTERING KERNELS AT CUBE CORNERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino-nucleon down and iso-e scattering kernals
!  Loop over incoming energy
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  e_in             = unu(j,k)
  e_in_l           = unub(j,k)
  e_in_u           = unub(j,k+1)

!-----------------------------------------------------------------------
!  Weak magnetism corrections for neutrino and antineutrino neutron and
!   proton scattering.
!-----------------------------------------------------------------------

  CALL nc_weak_mag( e_in, xi_p_wm, xi_n_wm, xib_p_wm, xib_n_wm )

!-----------------------------------------------------------------------
!  Neutrino-nucleon down and iso-e scattering kernals
!  Loop over ouytgoing energy
!-----------------------------------------------------------------------

  DO kp = 1,k
    e_out          = unu(j,kp)
    e_out_l        = unub(j,kp)
    e_out_u        = unub(j,kp+1)

!-----------------------------------------------------------------------
!  Compute neutrino-neutron elastic scattering kernels
!-----------------------------------------------------------------------

    sm             = mn
    
    CALL sctncal( e_in, e_out, e_in_l, e_in_u, e_out_l, e_out_u, tmev, cmpn, &
& sm, sct0_n, sct1_n, cv_n, ca_n )

!-----------------------------------------------------------------------
!  Compute neutrino-proton elastic scattering kernels
!-----------------------------------------------------------------------

    sm             = mp
    CALL sctncal( e_in, e_out, e_in_l, e_in_u, e_out_l, e_out_u, tmev, cmpp, &
& sm, sct0_p, sct1_p, cv_p, ca_p )

!-----------------------------------------------------------------------
!  Load scattering kernal array
!-----------------------------------------------------------------------

    sctn0 (j,k,kp,ida,ita,iya,i_ray) = REAL( LOG10( DABS(sct0_n * xi_n_wm  + sct0_p * xi_p_wm)  + epsilon ) )
    sctn0 (j,kp,k,ida,ita,iya,i_ray) = REAL( LOG10( DABS(sct0_n * xi_n_wm  + sct0_p * xi_p_wm)  + epsilon ) )
    sctnb0(j,k,kp,ida,ita,iya,i_ray) = REAL( LOG10( DABS(sct0_n * xib_n_wm + sct0_p * xib_p_wm) + epsilon ) )
    sctnb0(j,kp,k,ida,ita,iya,i_ray) = REAL( LOG10( DABS(sct0_n * xib_n_wm + sct0_p * xib_p_wm) + epsilon ) )

  END DO ! kp = 1,k
END DO ! k = 1,nnugpmx

RETURN
END SUBROUTINE scatnrgn

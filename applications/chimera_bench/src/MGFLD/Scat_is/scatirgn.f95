SUBROUTINE scatirgn( j, ij_ray, ik_ray, idd, itt, iyy, ida, ita, iya )
!-----------------------------------------------------------------------
!
!    File:         scatirgn
!    Module:       scatirgn
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/00
!
!    Purpose:
!      To call the subroutines that compute the isoenergetic
!       scattering rates and store the rates in arrays for
!       interpolation.
!
!    Subprograms called:
!        scatical
!
!    Input arguments:
!
!  j                    : radial zone number
!  ij_ray               : index denoting the j-index of a specific radial ray
!  ik_ray               : index denoting the k-index of a specific radial ray
!  in                   : 0, neutrino-neutron isoenergetic scattering omitted
!                         1, neutrino-neutron isoenergetic scattering included
!  ip                   : 0, neutrino-proton isoenergetic scattering omitted
!                         1, neutrino-proton isoenergetic scattering included
!  ihe                  : 0, neutrino-helium isoenergetic scattering omitted
!                         1, neutrino-helium isoenergetic scattering included
!  iheavy               : 0, neutrino-heavy nucleus isoenergetic scattering
!                           omitted
!                         1, neutrino-heavy nucleus isoenergetic scattering included
!  iscat                : 0, all neutrino scattering processes omitted
!                         1, neutrino scattering processes not necessarily omitted
!  idd                  : density grid index
!  itt                  : temperature grid index
!  iyy                  : electron fraction grid index
!  ida                  : absorption-emission function table density index
!  ita                  : absorption-emission function table temperature index
!  iya                  : absorption-emission function table electron fraction index
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  nnugp(n)             : number of energy zones for neutrinos of type n
!  unu(j,k)             : energy of energy zone k in radial zone j (MeV)
!  dunu(j,k)            : energy width of energy zone k in radial zone j (MeV)
!  dgrid(idty(j,ij_ray,ik_ray))
!                       : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and
!                          ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  estble(i,j,ida,ita,iya,ij_ray,ik_ray) : equation of state table array
!
!    Output arguments (common):
!
!  cohsct               : neutrino-nucleus (nucleon) isoenergetic scattering
!                          kernal array
!  cohbsct              : antineutrino-nucleus (nucleon) isoenergetic scattering
!                          kernal array
!
!    Include files:
!  kind_module, numerical_module
!  eos_snc_x_module, nu_dist_module, nu_energy_grid_module, scat_i_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, one, epsilon

USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, estble, escnst, idty
USE nu_dist_module, ONLY : unu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE scat_i_module, ONLY : cohsct, cohbsct

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! index denoting the k-index of a specific radial ray
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

REAL(KIND=double)                :: rho           ! density (g/cm**3)
REAL(KIND=double)                :: t             ! temperature (K)
REAL(KIND=double)                :: enu           ! neutrino zone-centered energy

REAL(KIND=double)                :: xn            ! neutron mass fraction
REAL(KIND=double)                :: xp            ! proton mass fraction
REAL(KIND=double)                :: xhe           ! helium mass fraction
REAL(KIND=double)                :: xh            ! heavy nucleus mass fraction
REAL(KIND=double)                :: ah            ! heavy nucleus mass number
REAL(KIND=double)                :: zh            ! heavy nucleus charge number

REAL(KIND=double)                :: rmdnns        ! mfp^-1 for neutrino-neutron scattering
REAL(KIND=double)                :: rmdnps        ! mfp^-1 for neutrino-proton scattering
REAL(KIND=double)                :: rmdnbns       ! mfp^-1 for antineutrino-neutron scattering
REAL(KIND=double)                :: rmdnbps       ! mfp^-1 for antineutrino-proton scattering
REAL(KIND=double)                :: rmdnhes       ! mfp^-1 for neutrino-helium scattering
REAL(KIND=double)                :: rmdnhs        ! mfp^-1 for neutrino-heavy nucleus scattering
REAL(KIND=double)                :: coh           ! sum of above for neutrinos
REAL(KIND=double)                :: cohb          ! sum of above for antineutrinos

REAL(KIND=double)                :: f10           ! function 10**

EXTERNAL f10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute independent variables at grid point
!-----------------------------------------------------------------------

rho                = f10   ( DBLE(idd)/dgrid(idty(j,ij_ray,ik_ray)) )
t                  = f10   ( DBLE(itt)/tgrid(idty(j,ij_ray,ik_ray)) )

!-----------------------------------------------------------------------
!  Fetch equation of state quantities needed for the computation of
!   the neutrino scattering functions.
!-----------------------------------------------------------------------

xn                 = f10( dble( estble(7 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 7,j,ij_ray,ik_ray) - epsilon
xp                 = f10( dble( estble(8 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 8,j,ij_ray,ik_ray) - epsilon
xh                 = f10( dble( estble(9 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst( 9,j,ij_ray,ik_ray) - epsilon
ah                 = f10( dble( estble(10,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst(10,j,ij_ray,ik_ray) - epsilon
zh                 = f10( dble( estble(11,j,ida,ita,iya,ij_ray,ik_ray) ) ) &
&                  - escnst(11,j,ij_ray,ik_ray) - epsilon

xn                 = DMAX1( xn, zero )
xp                 = DMAX1( xp, zero )
xh                 = DMAX1( xh, zero )
ah                 = DMAX1( ah, zero )
zh                 = DMAX1( zh, zero )
xhe                = DMAX1( one - xn - xp - xh, zero )

!-----------------------------------------------------------------------
!  Neutrino-nucleus (nucleon) isoenergetic scattering kernals
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  enu              = unu(j,k)
  CALL scatical( rho, t, enu, xn, xp, xhe, xh, ah, zh, rmdnns, rmdnps, &
&  rmdnbns, rmdnbps, rmdnhes, rmdnhs, coh, cohb )
  cohsct (j,k,ida,ita,iya,ij_ray,ik_ray) = REAL( DLOG10( coh  + epsilon ) )
  cohbsct(j,k,ida,ita,iya,ij_ray,ik_ray) = REAL( DLOG10( cohb + epsilon ) )
END DO

RETURN
END SUBROUTINE scatirgn

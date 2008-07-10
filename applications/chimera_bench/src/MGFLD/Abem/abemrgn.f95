SUBROUTINE abemrgn( j, ij_ray, ik_ray, idd, itt, iyy, ida, ita, iya )
!-----------------------------------------------------------------------
!
!    File:         abemrgn
!    Module:       abemrgn
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!                  7/19/07 WRH added tabular electron capture
!
!    Purpose:
!      To call the subroutines that compute the absorption and
!       emission inverse mean free paths and to store the total
!       inverse mean free paths in arrays for interpolation by
!       by subroutine abemrate.
!
!    Variables that must be passed through common:
!  nnugp(n)             : number of energy zones for neutrinos of type n
!  unu(k)               : energy of energy zone k (MeV)
!  dunu(k)              : energy width of energy zone k (MeV)
! dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
! tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
! ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)       : index for dgrid, tgrid, and ygrid for zone j
! estble(i,j,ida,ita,iya,ij_ray,ik_ray) : equation of state table array
!
!    Subprograms called:
!  abem_cal             : computes neutrino emission and absorption on free
!                          neutrons and protons
!  abemnc               : computes neutrino emission and absorption on heavy
!                          nuclei using the FFN formalizm
!  abemhnc              : computes neutrino emission and absorption on heavy
!                          nuclei using Haxton's rates
!
!    Input arguments:
!
!  j                    : radial zone number
!  ij_ray               : j-index of a radial ray
!  ik_ray               : k-index of a radial ray
!  idd                  : density grid index
!  itt                  : temperature grid index
!  iyy                  : electron fraction grid index
!  ida                  : absorption-emission function table density index
!  ita                  : absorption-emission function table temperature index
! iya                   : absorption-emission function table electron fraction index
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  iaefnp               : 0, e-neutrino and e-antineutrino absorption on free
!                          neutrons and protons omitted
!                       : 1, e-neutrino and e-antineutrino absorption on free 
!                          neutrons and protons included
!  rhoaefnp             : density above which e-neutrino and e-antineutrino
!                          absorption on free neutrons and protons is turned off
!  iaencnu              : 0, e-neutrino and e-antineutrino absorption on nuclei
!                          (Haxton's rates omitted
!                       : 1, e-neutrino and e-antineutrino absorption on 
!                          nuclei (Haxton's rates included
!  roaencnu             : density above which e-neutrino and e-antineutrino
!                          absorption on nuclei (Haxton's rates) is turned off
!  iaence               : 0, e-neutrino absorption on nuclei (FFN rates) omitted
!                         1, e-neutrino absorption on nuclei (FFN rates) included
!  edmpe                : ( en - cmpn ) - ( ep - cmpp) ; mean difference between
!                          the exitation energy of daughter and parent nucleus 
!                          for e-neutrino absorption on nuclei
!  iaenca               : 0, e-antineutrino absorption on nuclei (FFN rates) omitted
!                         1, e-antineutrinoneutrino absorption on nuclei (FFN rates) included
!  edmpa                : ( en - cmpn ) - ( ep - cmpp) ; mean difference between
!                          the exitation energy of daughter and parent nucleus
!                          for e-antineutrino absorption on nuclei
!  nnugp(n)             : number of energy zones for neutrinos of type n
!  unu(j,k)             : energy of energy zone k (MeV) in radial zone j
!  dgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in rho for zone j
!  tgrid(idty(j,ij_ray,ik_ray)) : number of table entries per decade in t for zone j
!  ygrid(idty(j,ij_ray,ik_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone j
!  idty(j,ij_ray,ik_ray)        : index for dgrid, tgrid, and ygrid for zone j
!  estble(i,j,ida,ita,iya,ij_ray,ik_ray)
!                       : equation of state table array
!  ye(j)                : electron fraction of zone j
!
!    Output arguments (common):
!
!   em(j,k,n,ida,ita,iya,ij_ray,ik_ray)
!                       : emission rate table array
!   ab(j,k,n,ida,ita,iya,ij_ray,ik_ray)
!                       : absorption rate table array
!
!    Modules used:
!  kind_module, array_module, numerical_module,
!  abem_module, edit_module, eos_snc_x_module, nu_dist_module,
!  nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu, ij_ray_dim
USE numerical_module, ONLY : zero, one, epsilon

USE abem_module, ONLY : em, ab
USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY : dgrid, tgrid, ygrid, estble, escnst, idty
USE nu_dist_module, ONLY : unu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j              ! shifted radial zone index
INTEGER, INTENT(in)              :: ij_ray         ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray         ! k-index of a radial ray
INTEGER, INTENT(in)              :: idd            ! density do index
INTEGER, INTENT(in)              :: itt            ! temperature do index
INTEGER, INTENT(in)              :: iyy            ! electron fraction do index
INTEGER, INTENT(in)              :: ida            ! density cube index
INTEGER, INTENT(in)              :: ita            ! temperature cube index
INTEGER, INTENT(in)              :: iya            ! electron fraction cube index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: k              ! incoming neutrino energy zone index
INTEGER                          :: n              ! neutrino flavor index

INTEGER                          :: i_ray          ! f(ij_ray,ik_ray)

REAL(KIND=double)                :: rho            ! density (g/cm**3)
REAL(KIND=double)                :: t              ! temperature (K)
REAL(KIND=double)                :: ye             ! electron fraction

REAL(KIND=double)                :: cmpn           ! neutron chemical potential
REAL(KIND=double)                :: cmpp           ! proton chemical potential
REAL(KIND=double)                :: cmpe           ! electron chemical potential
REAL(KIND=double)                :: xneut          ! free neutron mass fraction
REAL(KIND=double)                :: xprot          ! free proton mass fraction
REAL(KIND=double)                :: xh             ! heavy nuclei mass fraction
REAL(KIND=double)                :: ah             ! heavy nuclei mass number
REAL(KIND=double)                :: zh             ! heavy nuclei charge number
REAL(KIND=double)                :: w              ! neutrino energy

REAL(KIND=double)                :: emisnp         ! n-neutrino emission inv mfp on free nucleons
REAL(KIND=double)                :: absrnp         ! n-neutrino absorption inv mfp on free nucleons
REAL(KIND=double)                :: emisnc         ! n-neutrino emission inv mfp on nuclei (FFN)
REAL(KIND=double)                :: absrnc         ! n-neutrino absorption inv mfp on nuclei (FFN)
REAL(KIND=double)                :: emish          ! n-neutrino emission inv mfp on nuclei (Haxton)
REAL(KIND=double)                :: absrh          ! n-neutrino absorption inv mfp on nuclei (Haxton)
REAL(KIND=double)                :: emist(nnugpmx) ! n-neutrino emission inv mfp on nuclei (Table)
REAL(KIND=double)                :: absrt(nnugpmx) ! n-neutrino absorption inv mfp on nuclei (Table)


REAL(KIND=double)                :: f10            ! function 10**
EXTERNAL f10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute independent variables at grid point
!-----------------------------------------------------------------------

rho                = f10   ( DBLE(idd)/dgrid(idty(j,ij_ray,ik_ray)) )
t                  = f10   ( DBLE(itt)/tgrid(idty(j,ij_ray,ik_ray)) )
ye                 = one - ( DBLE(iyy)/ygrid(idty(j,ij_ray,ik_ray)) )

!-----------------------------------------------------------------------
!  Fetch equation of state quantities at the grid points needed for the
!   computation of the absorption and emission rates there.
!-----------------------------------------------------------------------

cmpn               = f10( DBLE( estble(4 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) - escnst( 4,j,ij_ray,ik_ray) - epsilon
cmpp               = f10( DBLE( estble(5 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) - escnst( 5,j,ij_ray,ik_ray) - epsilon
cmpe               = f10( DBLE( estble(6 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) - escnst( 6,j,ij_ray,ik_ray) - epsilon
xneut              = f10( DBLE( estble(7 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) - escnst( 7,j,ij_ray,ik_ray) - epsilon
xprot              = f10( DBLE( estble(8 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) - escnst( 8,j,ij_ray,ik_ray) - epsilon
xh                 = f10( DBLE( estble(9 ,j,ida,ita,iya,ij_ray,ik_ray) ) ) - escnst( 9,j,ij_ray,ik_ray) - epsilon
ah                 = f10( DBLE( estble(10,j,ida,ita,iya,ij_ray,ik_ray) ) ) - escnst(10,j,ij_ray,ik_ray) - epsilon
zh                 = f10( DBLE( estble(11,j,ida,ita,iya,ij_ray,ik_ray) ) ) - escnst(11,j,ij_ray,ik_ray) - epsilon

xneut              = DMAX1( xneut, zero )
xprot              = DMAX1( xprot, zero )
xh                 = DMAX1( xh   , zero )
ah                 = DMAX1( ah   , zero )
zh                 = DMAX1( zh   , zero )

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!
!          \\\\\ NEUTRINO ABSORPTION AND EMISSION RATES /////
!
!-----------------------------------------------------------------------

DO n = 1,nnu

!-----------------------------------------------------------------------
!  If no neutrinos, cycle
!-----------------------------------------------------------------------

  IF ( nnugp(n) == 0 ) CYCLE
!  IF ( n == 1  .and.  j == 2 ) WRITE(nlog,"(a,2i4,7es9.2)") &
!&  'Abem', j, n, t, ye, rho, xh, ah, xneut, xprot

!-----------------------------------------------------------------------
!  n-neutrino - nuclei absorption and emission inverse mean free
!   paths (/cm), from table for all energies.
!-----------------------------------------------------------------------

  CALL abemtnc( j, n, rho, t, ye, xh, ah, cmpn, cmpp, cmpe, absrt, emist )

  DO k = 1,nnugp(n)

    w              = unu(j,k)

!-----------------------------------------------------------------------
!  n-neutrino - free nucleon absorption and emission inverse mean
!   free paths (/cm)
!-----------------------------------------------------------------------

    CALL abem_cal( n, w, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
&    cmpe, absrnp, emisnp, ye )

!-----------------------------------------------------------------------
!  n-neutrino - nuclei absorption and emission inverse mean free
!   paths (/cm).
!-----------------------------------------------------------------------

    CALL abemnc( n, w, rho, t, xh, ah, zh, cmpn, cmpp, cmpe, absrnc, emisnc )

!-----------------------------------------------------------------------
!  n-neutrino - nuclei absorption and emission inverse mean free paths
!   (/cm), Haxton's rates.
!-----------------------------------------------------------------------

    CALL abemhnc( n, k, rho, t, xh, ah, cmpn, cmpp, cmpe, absrh, emish )

!-----------------------------------------------------------------------
!  Store total neutrino absorption and emission rates
!-----------------------------------------------------------------------

    em(j,k,n,ida,ita,iya,i_ray) = REAL( DLOG10( emisnp + emisnc + emish + emist(k) + epsilon ) )
    ab(j,k,n,ida,ita,iya,i_ray) = REAL( DLOG10( absrnp + absrnc + absrh + absrt(k) + epsilon ) )

!    If( n == 1  .and.  j == 2 ) WRITE(nlog,"(i2,9es10.3)") &
!&    k , w, emisnp, emisnc, emish, emist(k), absrnp, absrnc, absrh, absrt(k)

!........End k do loop

   END DO ! k = 1,nnugp(n)

!........End n do loop

END DO ! n = 1,nnu

RETURN
END SUBROUTINE abemrgn

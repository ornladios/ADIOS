SUBROUTINE abemrgn_z( k, ki_ray, kj_ray, idd, itt, iyy, ida, ita, iya, agr )
!-----------------------------------------------------------------------
!
!    File:         abemrgn_z
!    Module:       abemrgn_z
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To call the subroutines that compute the absorption and
!       emission inverse mean free paths and to store the total
!       inverse mean free paths in arrays for interpolation by
!       by subroutine abemrate.
!
!    Variables that must be passed through common:
!  nnugp(n)             : number of energy zones for neutrinos of type n
!  unui(ie)              : energy of energy zone ie (MeV)
!  dunu(ie)              : energy width of energy zone ie (MeV)
!  dgrid(idty(k,kj_ray,ki_ray)) : number of table entries per decade in rho for zone k
!  tgrid(idty(k,kj_ray,ki_ray)) : number of table entries per decade in t for zone k
!  ygrid(idty(k,kj_ray,ki_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone k
!  idty(k,kj_ray,ki_ray)        : index for dgrid, tgrid, and ygrid for zone k
!  estble(i,k,ida,ita,iya,kj_ray,ki_ray) : equation of state table array
!
!    Subprograms called:
!  abem_cal_z           : computes neutrino emission and absorption on free
!                          neutrons and protons
!
!    Input arguments:
!
!  k                    : azimuthal zone number
!  ki_ray               : x (radial) index of a specific z (azimuthal) ray
!  kj_ray               : y (angular) index of a specific z (azimuthal) ray
!  idd                  : density grid index
!  itt                  : temperature grid index
!  iyy                  : electron fraction grid index
!  ida                  : absorption-emission function table density index
!  ita                  : absorption-emission function table temperature index
!  iya                  : absorption-emission function table electron fraction index
!  agr                  : lapse function
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
!                       1, e-neutrino absorption on nuclei (FFN rates) included
!  edmpe                : ( en - cmpn ) - ( ep - cmpp) ; mean difference between
!                          the exitation energy of daughter and parent nucleus 
!                          for e-neutrino absorption on nuclei
!  iaenca               : 0, e-antineutrino absorption on nuclei (FFN rates) omitted
!                         1, e-antineutrinoneutrino absorption on nuclei (FFN rates) included
!  edmpa                : ( en - cmpn ) - ( ep - cmpp) ; mean difference between
!                          the exitation energy of daughter and parent nucleus
!                          for e-antineutrino absorption on nuclei
!  nnugp(n)             : number of energy zones for neutrinos of type n
!  unui(k,ie)             : energy of energy zone ie (MeV) in radial zone k
!  dgrid(idty(k,kj_ray,ki_ray)) : number of table entries per decade in rho for zone k
!  tgrid(idty(k,kj_ray,ki_ray)) : number of table entries per decade in t for zone k
!  ygrid(idty(k,kj_ray,ki_ray)) : number of table entries in ye between ye = 0.5 and ye = 0 for zone k
!  idty(k,kj_ray,ki_ray)        : index for dgrid, tgrid, and ygrid for zone k
!  estble(i,k,ida,ita,iya,kj_ray,ki_ray)
!                       : equation of state table array
!  ye(k)                : electron fraction of zone k
!
!    Output arguments (common):
!
!  em(k,ie,n,ida,ita,iya,kj_ray,ki_ray)
!                       : emission rate table array
!  ab(k,ie,n,ida,ita,iya,kj_ray,ki_ray)
!                       : absorption rate table array
!
!    Modules used:
!  kind_module, array_module, numerical_module,
!  abem_z_module, eos_snc_z_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nnu, ik_ray_dim
USE numerical_module, ONLY : zero, one, epsilon

USE abem_z_module, ONLY : em, ab
USE eos_snc_z_module, ONLY : dgrid, tgrid, ygrid, estble, escnst, idty
USE nu_energy_grid_module, ONLY : nnugp, unui

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: k             ! radial zone index
INTEGER, INTENT(in)              :: kj_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ki_ray        ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: idd           ! density do index
INTEGER, INTENT(in)              :: itt           ! temperature do index
INTEGER, INTENT(in)              :: iyy           ! electron fraction do index
INTEGER, INTENT(in)              :: ida           ! density cube index
INTEGER, INTENT(in)              :: ita           ! temperature cube index
INTEGER, INTENT(in)              :: iya           ! electron fraction cube index

REAL(KIND=double), INTENT(in)    :: agr           ! lapse function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ie            ! incoming neutrino energy zone index
INTEGER                          :: n             ! neutrino flavor index

INTEGER                          :: j_ray         ! f(kj_ray,ki_ray)

REAL(KIND=double)                :: rho           ! density (g/cm**3)
REAL(KIND=double)                :: t             ! temperature (K)
REAL(KIND=double)                :: ye            ! electron fraction

REAL(KIND=double)                :: cmpn          ! neutron chemical potential
REAL(KIND=double)                :: cmpp          ! proton chemical potential
REAL(KIND=double)                :: cmpe          ! electron chemical potential
REAL(KIND=double)                :: xneut         ! free neutron mass fraction
REAL(KIND=double)                :: xprot         ! free proton mass fraction
REAL(KIND=double)                :: xh            ! heavy nuclei mass fraction
REAL(KIND=double)                :: ah            ! heavy nuclei mass number
REAL(KIND=double)                :: zh            ! heavy nuclei charge number
REAL(KIND=double)                :: w             ! neutrino energy

REAL(KIND=double)                :: emisnp        ! n-neutrino emission inv mfp on free nucleons
REAL(KIND=double)                :: absrnp        ! n-neutrino absorption inv mfp on free nucleons
REAL(KIND=double)                :: emisnc        ! n-neutrino emission inv mfp on nuclei (FFN)
REAL(KIND=double)                :: absrnc        ! n-neutrino absorption inv mfp on nuclei (FFN)
REAL(KIND=double)                :: emish         ! n-neutrino emission inv mfp on nuclei (Haxton)
REAL(KIND=double)                :: absrh         ! n-neutrino absorption inv mfp on nuclei (Haxton)


REAL(KIND=double)                :: f10           ! function 10**
EXTERNAL f10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute independent variables at grid point
!-----------------------------------------------------------------------

rho                = f10   ( DBLE(idd)/dgrid(idty(k,kj_ray,ki_ray)) )
t                  = f10   ( DBLE(itt)/tgrid(idty(k,kj_ray,ki_ray)) )
ye                 = one - ( DBLE(iyy)/ygrid(idty(k,kj_ray,ki_ray)) )

!-----------------------------------------------------------------------
!  Fetch equation of state quantities at the grid points needed for the
!   computation of the absorption and emission rates there.
!-----------------------------------------------------------------------

cmpn               = f10( DBLE( estble(4 ,k,ida,ita,iya,kj_ray,ki_ray) ) ) - escnst( 4,k,kj_ray,ki_ray) - epsilon
cmpp               = f10( DBLE( estble(5 ,k,ida,ita,iya,kj_ray,ki_ray) ) ) - escnst( 5,k,kj_ray,ki_ray) - epsilon
cmpe               = f10( DBLE( estble(6 ,k,ida,ita,iya,kj_ray,ki_ray) ) ) - escnst( 6,k,kj_ray,ki_ray) - epsilon
xneut              = f10( DBLE( estble(7 ,k,ida,ita,iya,kj_ray,ki_ray) ) ) - escnst( 7,k,kj_ray,ki_ray) - epsilon
xprot              = f10( DBLE( estble(8 ,k,ida,ita,iya,kj_ray,ki_ray) ) ) - escnst( 8,k,kj_ray,ki_ray) - epsilon
xh                 = f10( DBLE( estble(9 ,k,ida,ita,iya,kj_ray,ki_ray) ) ) - escnst( 9,k,kj_ray,ki_ray) - epsilon
ah                 = f10( DBLE( estble(10,k,ida,ita,iya,kj_ray,ki_ray) ) ) - escnst(10,k,kj_ray,ki_ray) - epsilon
zh                 = f10( DBLE( estble(11,k,ida,ita,iya,kj_ray,ki_ray) ) ) - escnst(11,k,kj_ray,ki_ray) - epsilon

xneut              = DMAX1( xneut, zero )
xprot              = DMAX1( xprot, zero )
xh                 = DMAX1( xh   , zero )
ah                 = DMAX1( ah   , zero )
zh                 = DMAX1( zh   , zero )

j_ray              = ik_ray_dim * ( ki_ray - 1 ) + kj_ray

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

  DO ie = 1,nnugp(n)

    w              = unui(ie)/agr

!-----------------------------------------------------------------------
!  n-neutrino - free nucleon absorption and emission inverse mean free
!   paths (/cm)
!-----------------------------------------------------------------------

    CALL abem_cal_z( n, w, rho, t, xneut, xprot, xh, ah, zh, cmpn, cmpp, &
&    cmpe, absrnp, emisnp, ye )

!-----------------------------------------------------------------------
!  Store total neutrino absorption and emission rates
!-----------------------------------------------------------------------

    em(k,ie,n,ida,ita,iya,j_ray) = REAL( DLOG10( emisnp + emisnc + emish + epsilon ) )
    ab(k,ie,n,ida,ita,iya,j_ray) = REAL( DLOG10( absrnp + absrnc + absrh + epsilon ) )

!-----------------------------------------------------------------------
!  End ie do loop
!-----------------------------------------------------------------------

   END DO ! ie = 1,nnugp(n)

!-----------------------------------------------------------------------
!  End n do loop
!-----------------------------------------------------------------------

END DO ! n = 1,nnu

RETURN
END SUBROUTINE abemrgn_z

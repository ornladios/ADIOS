!-----------------------------------------------------------------------
!    Module:       e_advct_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE e_advct_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Neutrino energy advection controls
!-----------------------------------------------------------------------
!  ivc_x     : neutrino energy advection due to x-hydro switch.
!
!   ivc_x = 0 : x-hydro neutrino energy advection bypassed.
!   ivc_x = 1 : x-hydro neutrino energy advection included.
!
!  ivc_y     : neutrino energy advection due to y-hydro switch.
!
!   ivc_y = 0 : y-hydro neutrino energy advection bypassed.
!   ivc_y = 1 : y-hydro neutrino energy advection included.
!
!  ivc_z     : neutrino energy advection due to y-hydro switch.
!
!   ivc_z = 0 : z-hydro neutrino energy advection bypassed.
!   ivc_z = 1 : z-hydro neutrino energy advection included.
!
!  tau_advct : optical depth above which neutrinos are advected in
!   energy like a gamma = 4/3 gas
!
!  rhomin_y_eadvect : density below which y-neutrino energy advection
!   is turned off
!
!  rhomin_z_eadvect : density below which z-neutrino energy advection
!   is turned off
!-----------------------------------------------------------------------

INTEGER                                                    :: ivc_x
INTEGER                                                    :: ivc_y
INTEGER                                                    :: ivc_z

REAL(KIND=double)                                          :: tau_advct
REAL(KIND=double)                                          :: rhomin_y_eadvect
REAL(KIND=double)                                          :: rhomin_z_eadvect

!-----------------------------------------------------------------------
!  General edit parameters
!-----------------------------------------------------------------------
!  nprint : a parameter used in the call to an edit subroutine denoting 
!   the unit number of the file to which the print file is to be sent.
!-----------------------------------------------------------------------

INTEGER                                                    :: nprint

!-----------------------------------------------------------------------
!  Time step
!-----------------------------------------------------------------------
!  dtnph : the coordinate time step, set by hydro and nuclear reactions
!   (if the material is not in nse), between time cycle m and time
!   cycle m + 1.
!-----------------------------------------------------------------------

REAL(KIND=double)                                          :: dtnph

!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------
!  t_cntrl_e_advct(n): n-neutrino zero-moment change time step criterion
!   due to advection, i.e., maximum permitted
!     abs( dpsi0(j,k,n)/psi0(j,k,n) )(n=1,2,3).
!   due to advection.
!
!  dpsi0nmax(n) : maximum relative psi0 change arising from the energy
!   advection of neutrinos of type n.
!
!  jdpsi0nmax(n) : radiail zone limiting the time step due the relative
!   psi0 change arising from the energy advection of neutrinos of type n.
!
!  dpsivmx(n), psivmin(n) : parameters used in determining the psi0 change
!   time step due to compression, expansion and advection.
!
!     dpsivmx(n) = max( abs( dpsi0(j,k,n) )/( psi0(j,k,n) + psivmin(n) ) ) .
!
!  dt_nu_e_advct(n): maximum time step allowed by the neutrino energy
!   advection criterion for neutrinos of energy flavor n.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:)                         :: jdt_nu_e_advct

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: t_cntrl_e_advct
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dpsi0nmax
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: psivmin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dpsivmx
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dt_nu_e_advct

!-----------------------------------------------------------------------
!  Quantities needed to compute neutrino energy advection rate
!-----------------------------------------------------------------------
!  rjmh(j)   ; the zone-centered value of the circumferential radius at
!   time m.
!
!  rajmhj)   : the zone-centered value of the circumferential radius at
!   time m+1.
!
!  yjmh(j)   ; the zone-centered value of the y position of angular zone j
!   at time m.
!
!  yajmh(j)  : the zone-centered value of the y position of angular zone j
!   at time m+1.
!
!  zjmh(j)   ; the zone-centered value of the z position of azimuthal zone k
!   at time m.
!
!  zajmh(j)  : the zone-centered value of the z position of azimuthal zone k
!   at time m+1.
!
!  adot(j)   : the time derivative of agrjmh_(j) at time m+1/2.
!
!  adot_a(j) : adot(j)/( 0.5 * ( agrajmh(j) + agrjmh(j) ) ).
!
!  bdot_b(j) : the time derivative of bjmh(j) divided by bjmh(j) at
!   time m+1/2.
!
!  ddot(j)   : the time derivative of rho(j) at time m+1/2.
!
!  ddot_d(j) : the ddot(j)/( 0.5 * ( rhoi(j) + rho(j) ) ).
!
!  rdot(j)   ; the time derivative of rjmh(j) at time m+1/2.
!
!  rdot_r(j) : rdot(j)/( 0.5 * ( r(j) + ri(j) ) ).
!
!  v_e1(j)   : advection velocity of the radiation density.
!
!  v_e2(j)   : advection velocity of the radiation pressure.
!
!  v_e(j,k)  : total neutrino advection velocity for radial zone j-1/2,
!   energy zone k.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: rjmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: rajmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: yjmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: yajmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: zjmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: zajmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: adot
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: adot_a
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: bdot_b
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: ddot
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: ddot_d
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: rdot
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: rdot_r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: v_e1
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: v_e2
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: v_e

!-----------------------------------------------------------------------
!  More quantities needed to compute the neutrino energy advection
!-----------------------------------------------------------------------
!  nmin       : the number of first real zone (=1+6)
!
!  nmax       : the number of last  real zone (=nnugp(n)+6)
!
!  psi(k)     : neutrino occupation number : 1D array paddded with 6 ghost
!   zones.
!
!  xk(k)      : the energy grid zone edge locations.
!
!  dk(k)      :  xk(k+1) - xk(k).
!
!  xk_0(k)    : the energy grid zone edge locations prior to lagrangian
!   update.
!
!  dk_0(k)    : dk(k) prior to lagrangian update.
!
!  xk_1(k)    : the energy grid zone edge locations at time m+1.
!
!  dk_1(k)    : dk(k)  at time m+1.
!
!  vk(k)      : the total neutrino advection velocity: 1D array paddded
!   with 6 ghost zones.
!
!  dvolk(k)   : energy space volume.
!
!  dvolk_0(k) : energy space volume prior to lagrangian update.
!
!  dvolk_1(k) : energy space volume at time m+1.
!-----------------------------------------------------------------------

INTEGER                                                    :: nmin
INTEGER                                                    :: nmax

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: psi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: xk
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dk
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: xk_0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dk_0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: xk_1
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dk_1
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: vk
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dvolk
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dvolk_0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dvolk_1

!-----------------------------------------------------------------------
!  Quantities derived from the neutrino energy zoning
!-----------------------------------------------------------------------
!  ncoefa(j,k)  : the coefficient for computing the number density of
!   neutrinos at time m.
!
!  ncoefaa(j,k) : the coefficient for computing the number density of
!   neutrinos at time m+1.
!
!  ecoefa(j,k)  : the coefficient for computing the energy density of
!   neutrinos at the zone center at time m.
!
!  ecoefae(j,k) : the coefficient for computing the energy density of
!   neutrinos at the zone edge at time m.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: ncoefa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: ncoefaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: ecoefa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: ecoefaa

!-----------------------------------------------------------------------
!  State variables
!-----------------------------------------------------------------------
!  u(j)       : zone-centered velocity (cm s^{-1}). 
!
!  rho(j)     : the density of radial zone j at timestep m (g cm^{-3}). 
!
!  rhoa(j)    : the density of radial zone j at timestep m + 1 (g cm^{-3}). 
!
!  r(j)       : the radius of radial zone j at timestep m (cm).
!
!  ra(j)      : the radius of radial zone j at timestep m + 1 (cm).
!
!  y(j)       : the y posiiton of angular zone j at timestep m (cm).
!
!  ya(j)      : the y posiiton of angular zone j at timestep m + 1 (cm).
!
!  z(j)       : the z posiiton of azimuthal zone j at timestep m (cm).
!
!  za(j)      : the z posiiton of azimuthal zone j at timestep m + 1 (cm).
!
!  agrjmh(j)  : the lapse function ('00' metric component),
!     dt(proper) = agr*dt(coordinate)
!  defined at zone center j-1/2 at timestep m.
!
!  agrajmh(j) : the lapse function ('00' metric component) defined at zone
!   center j-1/2 at timestep m + 1, used in the neutrino transport modules.
!
!  rstmss(j)  : the rest mass enclosed by the outer boundary of radial
!   zone j (g).
!
!  dmrst(j)   : rstmss(j) - rstmss(j-1).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: u
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: rhoa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: ra
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: ya
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: za
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: agrjmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: agrajmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: rstmss
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: dmrst

!-----------------------------------------------------------------------
!  Neutrino diffusion
!-----------------------------------------------------------------------
!  psi0(j,k,n): the zero moment of the occupation distribution for neutrinos
!   at the midpoint of radial zone j, of energy zone k, and of type n.
!
!  psi1(j,k,n): the first moment of the occupation distribution for neutrinos
!   at the outer boundary radial zone j, of energy zone k, and of type n.
!
!  psi0_a(j,k,n): the zero moment of the occupation distribution for neutrinos
!   at the midpoint of radial zone j, of energy zone k, and of type m+1.
!
!  Ef(j,k,n): the Eddington factor at the outer boundary radial zone j,
!   of neutrino energy zone k, and of neutrino type n.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: psi0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: psi1
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: psi0_a
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: flxf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: Ef

!-----------------------------------------------------------------------
!  Energy advection parameters derived from psi0 and psi1
!-----------------------------------------------------------------------
!  unujv(j,n,ij_ray,ik_ray)    : the cumulative energy transferred to
!   n-type neutrinos by advection in radial zone j (ergs).
!
!  unucrv(n,ij_ray,ik_ray)     : the total energy transferred to n-type
!   neutrinos by advection (ergs).
!
!  dunujvdt(j,n,ij_ray,ik_ray) : the rate of energy transferred to
!   n-neutrinos by advection in radial zone j (ergs/s gram).
!
!  dndt_v(j,k,n,ij_ray,ik_ray) : net rate of n-neutrino production by
!   energy advection in zones j, k and neutrino type n (/cm3*sec*MeV).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: unucrv
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: dunujvdt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: unujv
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: dndt_v


END module e_advct_module

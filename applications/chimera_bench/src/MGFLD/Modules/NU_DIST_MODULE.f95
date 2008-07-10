!-----------------------------------------------------------------------
!    Module:       nu_dist_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE nu_dist_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Neutrino group energies
!-----------------------------------------------------------------------
!  unumn       : the minimum value of the zone-centered neutrino energy
!   (MeV). 
!
!  unumx       : the maximum value of the zone-centered neutrino energy
!   (MeV).
!
!  runu        : the ratio of adjecent energy zones (appropriate when
!   energy zoning is chosen to be geometric).
!     runu = unu(k+1)/unu(k)
!
!  unu(j,k)    : the value of the zone-centered neutrino energy of energy
!   zone k at radial zone center j at time m (MeV).
!
!  unub(j,k)   : the value of the inner zone-edge neutrino energy of
!   energy zone k at radial zone center j at time m (MeV).
!
!  dunu(j,k)   : the width of energy zone k at radial zone center j (MeV),
!   dunu(j,k) = unub(j,k+1) - unub(j,k) at time m
!
!  unue(j,k)   : the value of the zone-centered neutrino energy of energy
!   zone k at radial zone edge j at time m (MeV).
!
!  unube(j,k)  : the value of the inner zone-edge neutrino energy of
!   energy zone k at radial zone edge j at time m (MeV).
!
!  dunue(j,k)  : the width of energy zone k at radial zone edge j (MeV),
!   dunu(j,k) = unub(j,k+1) - unub(j,k)
!
!  stwt(n)     : the statistical weight of neutrinos of type n, i.e.,
!   the number of distinct multiplicities of neutrinos subsumed as
!   neutrinos of type n. (Typically, stwt(1) = 1., stwt(2) = 1.,
!   stwt(3) = 3. and stwt(3) = 2.)
!-----------------------------------------------------------------------

REAL(KIND=double)                                     :: unumn
REAL(KIND=double)                                     :: unumx
REAL(KIND=double)                                     :: runu

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: unu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: dunu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: unub

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: unue
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: dunue
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: unube

!-----------------------------------------------------------------------
!  Quantities derived from the neutrino energy zoning
!-----------------------------------------------------------------------
!  ncoefa(j,k)   : the coefficient for computing the number density of
!   neutrinos at time m.
!
!  ecoefa(j,k)   : the coefficient for computing the energy density of
!   neutrinos at the zone center at time m.
!
!  ecoefae(j,k)  : the coefficient for computing the energy density of
!   neutrinos at the zone edge at time m.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: ncoefa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: ecoefa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: ecoefae

!-----------------------------------------------------------------------
!  Neutrino occupation distribution
!-----------------------------------------------------------------------
!  psi0(j,k,n) : the zero moment of the occupation distribution for
!   neutrinos at the midpoint of radial zone j, of energy zone k, and
!   of type n.
!
!  psi1(j,k,n) : the first moment of the occupation distribution for
!   neutrinos at the outer boundary radial zone j, of energy zone k,
!   and of type n.
!
!  stwt(n)     : the number of neutrino flavors grouped as type n
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: psi0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: psi1

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: stwt

!-----------------------------------------------------------------------
!  Neutrino number density
!-----------------------------------------------------------------------
!  n_nu(j,k,n) : the number of neutrinos per unit volume at the midpoint
!   of radial zone j, of energy zone k, and of type n.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: n_nu

!-----------------------------------------------------------------------
!  Neutrino diffusion
!-----------------------------------------------------------------------
!  rhs1(j,k,n)      : the right-hand side of the first moment of the
!   transport equation divided by psi1(j,k,n)/psi0(j,k,n)
!
!  tmfp(j,k,n)      : the transport mean free path for neutrinos of
!   energy zone k, type n, at radial zone j-1/2 and at time m + 1.
!
!  tmfp_j(j,k,n)    : the transport mean free path for neutrinos of
!   energy zone k, type n, at radial zone j and at time m + 1.
!
!  e_mfp_inv(j,k,n) : the energy inverse mean free path for neutrinos of
!   energy zone k, type n, at radial zone j-1/2 and at time m + 1.
!
!  dc(j,k,n)        : the diffusion coefficient for neutrinos of energy
!   zone k, type n, at radial zone j and at time m + 1.
!
!  dcr(j,k,n)       : the diffusion coefficient for neutrinos of energy
!   zone k, type n, at radial zone j and at time m.
!
!  fluxlm(j,k,n)    : the flux limiter for neutrinos of energy zone k,
!   type n, at radial zone j and at time m + 1.
!
!  comvcf(j,k,n)    : the matter-neutrino decoupling parameter.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: rhs1
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: tmfp
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: tmfp_j
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: e_mfp_inv
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: dc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: dcr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: fluxlm

!-----------------------------------------------------------------------
!  State variables for use in the transport subroutines
!-----------------------------------------------------------------------
!  rjmh(j)       : the radius of radial zone j-1/2 at timestep m + 1 used
!   in the neutrino modules (g/cm**3).
!
!  gamgr_nu(j)   : the inverse radial metric component for zone j at
!   timestep m used in the neutrino transport modules.
!
!  agr_nu(j)     : the lapse function ('00' metric component),
!     dt(proper) = agr*dt(coordinate)
!   transport modules. defined at zone boundary j at timestep m, used in
!   the neutrino
!
!  agrjmh_nu(j)  : the lapse function ('00' metric component),
!     dt(proper) = agr*dt(coordinate)
!   defined at zone center j-1/2 at timestep m, used in the neutrino
!   transport modules.
!
!  area(j)       : the area of radial zone j (cm**2).
!
!  areajmh(j)    : the area of radial zone j-1/2 (cm**2).
!
!  vol(j)        :  c * dtngrjmh(j) proper volume of computational
!   volume between r(j) and r(j-1) (cm**3).
!
!  drjmh(j)      : the proper distance between r(j+1/2) and r(j-1/2) (cm).
!
!  drjmh_inv(j)  : 1/drjmh(j).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: rjmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: gamgr_nu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: agr_nu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: agrjmh_nu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: area
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: areajmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: vol
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: drjmh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: drjmh_inv

!-----------------------------------------------------------------------
!  Thermodynamic quantities derived from the neutrino occupation
!   distribution
!-----------------------------------------------------------------------
!  aunu(j)   : the neutrino energy/mass (ergs/g) of radial zone j.
!
!  apnu(j)   : the neutrino pressure (dynes/cm2) of radial zone j.
!
!  apnun(j)  : the neutrino scalar pressure in radial zone j at the
!   last neutrino time cycle (dynes/cm2).
!
!  apnur(j)  : the neutrino scalar pressure in radial zone j at the
!   last minus 1 neutrino time cycle (dynes/cm2).
!
!  apmnun(j) : the material pressure in radial zone j at the last
!   neutrino time cycle (dynes/cm2).
!
!  apmnur(j) : the material pressure in radial zone j at the last
!   minus 1 neutrino time cycle (dynes/cm2).
!
!  rhonun(j) : the material density in radial zone j at the last
!   neutrino time cycle (g/cm3).
!
!  rhonur(j) : the material density in radial zone j at the last
!   minus 1 neutrino time cycle (g/cm3).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: aunu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: apnu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: apnun
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: apnur
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: apmnun
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: apmnur
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: rhonun
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: rhonur

!-----------------------------------------------------------------------
!  Fluxes and Stresses
!-----------------------------------------------------------------------
!
!  fnu(j)          : the energy flux of all neutrinos across radial zone
!   boundary j (ergs/cm2*sec).
!
!  fluxnuk(j,k,n)  : the energy flux of n-neutrinos of group k
!   across radial zone j.
!
!  fluxnu(j,n)     : the energy flux of n-neutrinos across of radial
!   zone j.
!
!  strsnu_x(j,k,n) : the force per unit mass exerted on the matter at
!   radial zone interface j by n-neutrinos of energy zone k (dynes/g).
!
!  stress_x(j,n,ij_ray,ik_ray) : the force per unit mass exerted on the
!   matter at radial zone int erface j by all n-neutrinos (dynes/g).
!
!  nu_str_ex(j)    : total force per unit mass exerted on the matter
!   at radial zone interface j by all n-neutrinos (dynes/g).
!
!  nu_str_cx(j)    : total force per unit mass exerted on the matter
!   at radial zone center j by all n-neutrinos (dynes/g).
!
!  strsnu_y(j,k,n) : the force per unit mass exerted on the matter at
!   angular zone interface j by n-neutrinos of energy zone k (dynes/g).
!
!  stress_y(j,n,j_ray,ik_ray) : the force per unit mass exerted on the
!   matter at angular zone interface j by all n-neutrinos (dynes/g).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: fnu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: fluxnuk
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: fluxnu
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: strsnu_x
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: stress_x
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: nu_str_ex
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: nu_str_cx
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: strsnu_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: stress_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: strsnu_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: stress_z

!-----------------------------------------------------------------------
!  Nwutrinospheres
!-----------------------------------------------------------------------
!  j_sphere(k,n,ij_ray,ik_ray) : the radial zone of the neutrinosphere
!   for neutrinos of type n and energy k (cm).                                   
!
!  r_sphere(k,n,ij_ray,ik_ray) : the radius of the neutrinosphere for
!   neutrinos of type n and energy k (cm).
!
!  d_sphere(k,n,ij_ray,ik_ray) : the density of the neutrinosphere for
!   neutrinos of type n and energy k (cm).
!
!  t_sphere(k,n,ij_ray,ik_ray) : the temperature of the neutrinosphere
!   for neutrinos of type n and energy k (cm).
!
!  m_sphere(k,n,ij_ray,ik_ray) : the mass enclosed by the neutrinosphere
!   of neutrinos of type n and energy k (cm).
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:)              :: j_sphere

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: r_sphere
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: d_sphere
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: t_sphere
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: m_sphere

!-----------------------------------------------------------------------
!  Transport parameters derived from psi0 and psi1
!-----------------------------------------------------------------------
!  dnurad(j,k,n,ij_ray,ik_ray) : the total number of n-neutrinos of
!   energy zone k per unit energy that have crossed the outer boundary
!   radial zone j (/MeV).
!  
!  unulum(n)                   : the core n-neutrino luminosity
!   (ergs/sec).
!
!  e_rad(ij_ray,ik_ray)        : the cumulative material energy entering
!   (negative) or leaving (positive) the grid (ergs).
!
!  unurad(n,ij_ray,ik_ray)     : the cumulative energy emitted from the
!   core in n-type neutrinos (ergs).
!
!  unucr(n)                    : the total energy of n-type neutrinos
!   currently residing in the core (ergs).
!
!  unuinfty(n)                 : the total energy of n-type neutrinos
!   currently residing in the core as observed at rest at infinity (ergs).
!
!  unucrt(n,ij_ray,ik_ray)     : the total energy transferred to
!   n-neutrinos by microphysics (ergs).
!
!  unukrad(k,n,ij_ray,ik_ray)  : the cumulative energy emitted from the
!   core by n-type neutrinos of energy zone k (ergs).
!
!  unujrad(j,n,ij_ray,ik_ray)  : the cumulative energy transported across
!   the outer boundary of radial zone j by n-type neutrinos (ergs).
!
!  unujcr(j,n)                 : the total energy of n-type neutrinos
!   currently residing in radial zone j (ergs).
!
!  unujinfty(j,n)              : the total energy of n-type neutrinos
!   currently residing in radial zone j as observed at rest at infinity
!   (ergs).
!
!  nnulum(n)                   : the core n-neutrino number luminosity
!
!  elec_rad(ij_ray,ik_ray)     : the net number of electrons advected in
!   (negative) or out (positive) of the grid
!
!  nnurad(n,ij_ray,ik_ray)     : the cumulative number of n-type
!   neutrinos emitted by the core.
!
!  nnucr(n,ij_ray,ik_ray)      : the net number of n-type neutrinos
!   currently residing in the core.
!
!  nnukrad(k,n,ij_ray,ik_ray)  : the cumulative number of n-type
!   neutrinos of energy k emitted by the core.
!
!  nnujrad(j,n,ij_ray,ik_ray)  : the net number of n-type neutrinos
!   transported across the outer boundary of radial zone j.
! 
!  nnujcr(j,n)                 : the net number of n-type neutrinos
!   currently residing n radial zone j.
!
!  dunujtdt(j,n)               :  the total rate of energy transferred
!   to n-neutrinos by microphysics in radial zone j (ergs/s gram).
!
!  unujt(j,n,ij_ray,ik_ray)    : the total cumulative energy transferred
!  xv n-type neutrinos by microphysics in radial zone j (ergs).
!
!  dudt_nu(j,ij_ray,ik_ray)    : energy deposition rate by all neutrinos
!   in radial zone j (ergs g^{-1} s^{-1}).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)  :: dnurad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: unulum
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: e_rad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unurad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: unucr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: unuinfty
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unucrt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unukrad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujrad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: unujcr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: unujinfty
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)          :: nnulum
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: elec_rad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: nnurad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: nnucr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: nnukrad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: nnujrad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: nnujcr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: dunujtdt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: dudt_nu

!-----------------------------------------------------------------------
!  Absorption and emission parameters derived from psi0 and psi1
!-----------------------------------------------------------------------
!  unucrea(nj_ray)              : the total energy transferred to n-type
!   neutrinos by emission and absorption (ergs).
!
!  dunujeadt(j,n,ij_ray,ik_ray) : the rate of energy transferred to n
!   neutrinos by emission and absorption in radial zone j (ergs/s gram).
!
!  unujea(j,n,ij_ray,ik_ray)    : the cumulative energy transferred to
!   n-type neutrinos by emission and absorption in  radial zone j (ergs).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unucrea
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: dunujeadt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujea

!-----------------------------------------------------------------------
!  Neutrino-electron scattering parameters derived from psi0 and psi1
!-----------------------------------------------------------------------
!  unucrnis(n,ij_ray,ik_ray)     : the total energy transferred to
!   n-type neutrinos by neutrino-electron scattering (ergs).
!
!  dunujnisdt(j,n,ij_ray,ik_ray) : the rate of energy transferred to
!   n-neutrinos by neutrino-electron scattering in radial zone j
!   (ergs/s gram).
!
!  unujnis(j,n,ij_ray,ik_ray)    : the cumulative energy transferred to
!   n-type neutrinos by neutrino-electron scattering in radial zone j
!   (ergs).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unucrnis
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: dunujnisdt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujnis

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation parameters derived from psi0 and
!   psi1
!-----------------------------------------------------------------------
!  unucrpa(n,ij_ray,ik_ray)     : the total energy transferred to n-type
!   neutrinos by electron-positron pair-annihilation (ergs).
!
!  dunujpadt(j,n,ij_ray,ik_ray) : the rate of energy transferred to
!   n-neutrinos by electron-positron pair-annihilation in radial zone j
!   (ergs/s gram).
!
!  unujpa(j,n,ij_ray,ik_ray)    : the cumulative energy transferred to
!   n-type neutrinos by electron-positron pair-annihilation in radial
!   zone j (ergs).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unucrpa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: dunujpadt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujpa

!-----------------------------------------------------------------------
!  Bremsstrahlung parameters derived from psi0 and psi1
!-----------------------------------------------------------------------
!  unucrba(n,ij_ray,ik_ray)     : the total energy transferred to n-type
!   neutrinos by bremsstrahlung pair-annihilation (ergs).
!
!  dunujbadt(j,n,ij_ray,ik_ray) : the rate of energy transferred to
!   n-neutrinos by bremsstrahlung pair-annihilation in radial zone j
!  (ergs/s gram).
!
!  unujba(j,n,ij_ray,ik_ray)    : the cumulative energy transferred to
!   n-type neutrinos by bremsstruhlung pair-annihilationin radial zone j
!   (ergs).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unucrba
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: dunujbadt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujba

!-----------------------------------------------------------------------
!  Neutrino-nucleus elastic scattering parameters derived from psi0 and
!   psi1
!-----------------------------------------------------------------------
!  unucrns(n,ij_ray,ik_ray)     : the total energy transferred to n-type
!   neutrinos by neutrino-nucleon inelastic scattering (ergs).
!
!  dunujnsdt(j,n,ij_ray,ik_ray) : the rate of energy transferred to
!   n-neutrinos by neutrino-nucleon inelastic scattering in radial zone j
!   (ergs/s gram).
!
!  unujns(j,n,ij_ray,ik_ray)    : the cumulative energy transferred to
!   n-type neutrinos by neutrino-nucleon inelastic scattering in radial
!   zone j (ergs).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unucrns
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: dunujnsdt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujns

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering parameters derived from psi0
!   and psi1
!-----------------------------------------------------------------------
!  unucrnns(n,ij_ray,ik_ray)     : the total energy transferred to n-type
!   neutrinos by neutrino-nucleon inelastic scattering (ergs).
!
!  dunujnnsdt(j,n,ij_ray,ik_ray) : the rate of energy transferred to
!   n-neutrinos by neutrino-nucleon inelastic scattering in radial zone j
!   (ergs/s gram).
!
!  unujnns(j,n,ij_ray,ik_ray)    : the cumulative energy transferred to
!   n-type neutrinos by neutrino-nucleon inelastic scattering in radial
!   zone j (ergs).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unucrnns
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: dunujnnsdt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujnns

!-----------------------------------------------------------------------
!  Parameterization of the neutrino distribution
!-----------------------------------------------------------------------
!  enuta(j,n)  : the temperature of a Fermi-Dirac neutrino distribution
!   having the same first three moments as the n-neutrino distribution
!   (MeV).
!
!  enucpa(j,n) : the chemical potential of a Fermi-Dirac neutrino
!   distribution having the same first three moments as the n-neutrino
!   distribution (MeV).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: enuta
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: enucpa

!-----------------------------------------------------------------------
!  Composition mass fractions
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: xn


END module nu_dist_module

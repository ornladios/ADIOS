!-----------------------------------------------------------------------
!    Module:       mdl_cnfg_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE mdl_cnfg_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Radial mass zoning
!-----------------------------------------------------------------------
!  jr_min : the outermost radial zone.
!  jr_max : the outermost radial zone.
!  jmr: the outermost radial zone of initial model just prior to rezoning.
!-----------------------------------------------------------------------

INTEGER                                        :: jr_min
INTEGER                                        :: jr_max
INTEGER                                        :: jmr

!-----------------------------------------------------------------------
!  Thermodynamic state
!-----------------------------------------------------------------------
!  rho(j)   : the density of radial zone j at timestep m (current time)
!   (g cm{-3}).
!
!  t(j)     : the temperature of radial zone j at timestep m (K).
!
!  ye(j)    : the electron fraction of radial zone j at timestep m.
!
!  rhor(j)  : the density of radial zone j at timestep m - 1 (g/cm**3).
!
!  tr(j)    : the temperature of radial zone j at timestep m - 1 (K).
!
!  yer(j)   : the electron fraction of radial zone j at timestep m - 1.
!
!  rho_i(j) : the initial density of radial zone j at timestep m (g cm{-3}).
!
!  t_i(j)   : the initial temperature of radial zone j at timestep m (K).
!
!  ye_i(j)  : the initial electron fraction of radial zone j at timestep m.
!
!  rho0(j)  : the density of radial zone j at problem initiation (g cm{-3}).
!
!  ye0(j)   : the electron fraction of radial zone j at problem initiation.
!
!  p0(j)    : the pressure of radial zone j at problem initiation
!   (dynes cm{-2}).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: t
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: ye
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rhor
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: tr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: yer
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rho0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: ye0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: p0

!-----------------------------------------------------------------------
!  Mechanical state
!-----------------------------------------------------------------------
!  u(j)    : the zone-centered x-velocity velocity radial zone j at
!   timestep m (cm s^{-1}).
!
!  v(j)    : the zone-centered y-velocity velocity radial zone j at
!   timestep m (cm s^{-1}).
!
!  w(j)    : the zone-centered z-velocity velocity radial zone j at
!   timestep m (cm s^{-1}).
!
!  r(j)    : the radius of outer boundary of radial zone j at timestep m
!   (cm).
!
!  dr(j) = r(j) - r(j-1).
!
!  r_i(j)  : the initial radius of outer boundary of radial zone j at
!   timestep m + 1 (cm).
!
!  dr_i(j) : r_i(j) - r_i(j-1).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: u
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: v
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: w
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dr

!-----------------------------------------------------------------------
!  Zone masses
!-----------------------------------------------------------------------
!  rstmss(j) : the rest mass enclosed by the outer boundary of radial
!   zone j (g).
!
!  dmrst(j)  : rstmss(j) - rstmss(j-1).
!
!  grvmss(j) : the gravitational mass enclosed by the outer boundary of
!   radial zone j at timestep
!       m (g).
!
!  dmgrv(j)  : grvmss(j) - grvmss(j-1).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rstmss
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dmrst
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: grvmss
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dmgrv

!-----------------------------------------------------------------------
!  General relativistic parameters
!-----------------------------------------------------------------------
!  gamgr(j)    : the inverse radial metric component,
!     gamgr * dr(proper) = dr(coordinate)
!   defined at zone boundary j at timestep m.
!
!  agrr(j)     : the lapse function or '00' metric component,
!     dt(proper) = agr * dt(coordinate),
!   defined at zone boundary j at timestep m - 1.
!
!  agr(j)      : the lapse function or '00' metric component,
!     dt(proper) = agr * dt(coordinate),
!   defined at zone boundary j at timestep m.
!
!  agrh(j)     : the lapse function or '00' metric component,
!     dt(proper) = agr * dt(coordinate),
!   defined at zone center j-1/2 at timestep m.
!
!  agra(j)     : the lapse function or '00' metric component,
!     dt(proper) = agr * dt(coordinate),
!   defined at zone boundary j at timestep m after Lagrangian update.
!
!  agrah(j)    : the lapse function or '00' metric component,
!     dt(proper) = agr * dt(coordinate),
!   defined at zone center j-1/2 at timestep m after Lagrangian update.
!
!  wgrr(j)     : the relativistic enthalpy for zone j at timestep m - 1.
!
!  wgr(j)      : the relativistic enthalpy for zone j at timestep m.
!
!  gamgrr(j)   : the inverse radial metric component for zone j at
!   timestep m - 1.
!
!  gamgra(j)   : the inverse radial metric component for zone j at
!   timestep m + 1.
!
!  grvmssa(j)  : the gravitational mass enclosed by the outer boundary
!   of radial zone j at timestep
!   m + 1 (g).
!
!  dmgrva(j)   : grvmssa(j) - grvmssa(j-1).
!
!  gamgr_i(j)  : the initial inverse radial metric component for zone j
!   at timestep m.
!
!  agr_i(j)    : the intial lapse function or '00' metric component,
!     dt(proper) = agr * dt(coordinate),
!   defined at zone boundary j at timestep m.
!
!  wgr_i(j)    : the initial relativistic enthalpy for zone j at
!   timestep m.
!
!  grvmss_i(j) : the initial gravitational mass enclosed by the outer
!   boundary of radial zone 
!   j at timestep m + 1 (g).
!
!  dmgrv_i(j   : grvmss_i(j) - grvmss_i(j-1).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: gamgr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: agrr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: agr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: agrh
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: agra
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: agrah
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: wgrr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: wgr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: gamgrr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: gamgra
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: grvmssa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dmgrva
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: gamgr_i
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: agr_i
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: wgr_i
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: grvmss_i
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dmgrv_i

!-----------------------------------------------------------------------
!  Maximum abd minimum central densities
!-----------------------------------------------------------------------
!  rojmax : the maximum density attained by inner radial zone between
!   edits of model configuration (via calls to subroutine editc) (g/cm**3).
!
!  rojmin : the minimum density attained by inner radial zone between
!   edits of model configuration (via calls to subroutine editc) (g/cm**3).
!-----------------------------------------------------------------------

REAL(KIND=double)                              :: rojmax
REAL(KIND=double)                              :: rojmin

!-----------------------------------------------------------------------
!  Specific energy density arrays for rezoning
!-----------------------------------------------------------------------
!  uad(j)   : the specific energy added to radial zone j to set the GR
!   gravity equal to the Newtonian gravity at problem initiation (ergs/g).
!
!  utoti(j) : the total energy per unit mass, computed before rezoning
!   (ergs/g).
!
!  utotf(j) : the total energy per unit mass, computed after rezoning
!   (ergs/g).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: uad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: utoti
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: utotf


END module mdl_cnfg_module

SUBROUTINE nu_energy_advct_inout_x( imin, imax, nx, ij_ray, ik_ray, &
& ij_ray_dim, ik_ray_dim, nez, nnu, nprintp, rhoip, rhop, rip, rp, tp, &
& yep, up, psi0p, psi1p, dtime, agr_c, dtp, jdtp )
!-----------------------------------------------------------------------
!
!    File:         nu_energy_advct_inout_x
!    Module:       nu_energy_advct_inout_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To receive variables from radial_ray_module, perform the neutrino energy,
!         advection step, and return to radial_ray_module the updated variables.
!
!    Input arguments:
!
!  imin              : minimum x-array index
!  imax              : maximum x-array index
!  nx                : x-array extent
!  ij_ray            : index denoting the j-index of a specific radial ray
!  ik_ray            : index denoting the k-index of a specific radial ray
!  ij_ray_dim        : number of y-zones on a processor before swapping
!  ik_ray_dim        : number of z-zones on a processor before swapping
!  nez               : neutrino energy array extent
!  nnu               : neutrino flavor array extent
!  nprintp           : unit number to print diagnostics
!  rhoip             : density before hydro advance (cm^{-3})
!  rhop              : density after hydro advance (cm^{-3})
!  rip               : radial zone radii before hydro advance (cm)
!  rp                : radial zone radii after hydro advance (cm)
!  tp                : temperature (K)
!  yep               : electron fraction
!  up                : zone-centered x-velocity velocity (cm s^{-1})
!  psi0p             : initial zero angular moments of the neutrino occupation number
!  psi1p             : first angular moments of the neutrino occupation number
!  dtime             : time step
!  agr_c             : zone-centered value of the lapse function
!
!    Output arguments:
!
!  psi0p             : updated zero angular moments of the neutrino occupation number
!  psi1p             : first angular moments of the neutrino occupation number
!  dtp               : minimum e-advection time step restrictions (s)
!  jdtp              : zones causing minimum time step restrictions dtp
!
!    Subprograms called:
!  nu_energy_advct_x : performs the x neutrino energy advection
!  psi1_cal          : computes the psi1's based on the updated psi0's
!  e_advct_time_step : computes the neutrino energy advection time step restrictions
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith

USE e_advct_module, ONLY : ivc_x, rho, rhoa, r, ra, agrjmh, agrajmh, u, &
& dmrst, rstmss, psi0, psi1, unujv, unucrv, jdt_nu_e_advct, dt_nu_e_advct, &
& dpsi0nmax, dtnph, nprint
USE nu_energy_grid_module, ONLY : nnugpmx
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin             ! minimum x-array index
INTEGER, INTENT(in)              :: imax             ! maximum x-array index
INTEGER, INTENT(in)              :: nx               ! lx-array extent

INTEGER, INTENT(in)              :: ij_ray           ! iindex denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray           ! iindex denoting the k-index of a specific radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping

INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent

INTEGER, INTENT(in)              :: nprintp          ! unit number to print diagnostics

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: rhoip    ! density before Lagrangian hydro step (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: rhop     ! density after Lagrangian hydro step (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                       :: rip      ! radial zone radii before hydro advance (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: rp       ! radial zone radii after hydro advance (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: tp       ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: yep      ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: up       ! zone-centered x-velocity velocity (cm s^{-1})
REAL(KIND=double), INTENT(in)                                        :: dtime    ! time step
REAL(KIND=double), INTENT(in), DIMENSION(nx  ,ij_ray_dim,ik_ray_dim) :: agr_c    ! zone-centered value of the lapse function

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER,           INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)  :: jdtp ! zone causing dtp

REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)  :: dtp  ! minimum allowed time step

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0p ! zero moment of the neutrino occupation probability
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1p ! first moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j               ! radial zone index
INTEGER                          :: jr_min          ! minimum radial zone index
INTEGER                          :: jr_max          ! maximum radial zone index

REAL(KIND=double), DIMENSION(nx) :: t               ! temperature (K)
REAL(KIND=double), DIMENSION(nx) :: ye              ! electron fraction

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Edit unit number
!-----------------------------------------------------------------------

nprint                       = nprintp

!-----------------------------------------------------------------------
!
!             \\\\\ SET UP FOR NEURINO ENERGY ADVECTION ///
!
!  Load variables received from radial_ray_module into e_advct_module
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize zone-center array size
!-----------------------------------------------------------------------

jr_min                       = imin + 1
jr_max                       = imax + 1

!-----------------------------------------------------------------------
!  Initialize time step parameters
!-----------------------------------------------------------------------

dpsi0nmax                    = zero
jdt_nu_e_advct               = 0
dt_nu_e_advct                = 1.d+20

!-----------------------------------------------------------------------
!  Return if nnugpmx = 0  or  ivc_x = 0
!-----------------------------------------------------------------------

IF ( nnugpmx == 0  .or.  ivc_x == 0 ) THEN
  dtp (40+1:40+nnu,ij_ray,ik_ray) = dt_nu_e_advct (1:nnu)
  jdtp(40+1:40+nnu,ij_ray,ik_ray) = jdt_nu_e_advct(1:nnu)
  RETURN
END IF ! nnugpmx == 0  .or.  ivc_x == 0

!-----------------------------------------------------------------------
!  Transfer psi0 to e_advct_module
!-----------------------------------------------------------------------

psi0(jr_min:jr_max,:,:)      = psi0p(imin:imax,:,:,ij_ray,ik_ray)
psi1(imin:imax+1,:,:)        = psi1p(imin:imax+1,:,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer state variables to shifted arrays in e_advct_module
!  No change in lapse so laod both agrjmh and agrajmh from agr_c
!-----------------------------------------------------------------------

rho    (jr_min:jr_max)       = rhoip     (imin:imax,ij_ray,ik_ray)
rhoa   (jr_min:jr_max)       = rhop      (imin:imax,ij_ray,ik_ray)
t      (jr_min:jr_max)       = tp        (imin:imax,ij_ray,ik_ray)
ye     (jr_min:jr_max)       = yep       (imin:imax,ij_ray,ik_ray)
u      (jr_min:jr_max)       = up        (imin:imax,ij_ray,ik_ray)
agrjmh (jr_min:jr_max)       = agr_c     (imin:imax,ij_ray,ik_ray)
agrajmh(jr_min:jr_max)       = agr_c     (imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer x-coordinates to e_advct_module
!-----------------------------------------------------------------------

r (imin:imax+1)              = rip(imin:imax+1)
ra(imin:imax+1)              = rp (imin:imax+1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Newtonian rest masses
!-----------------------------------------------------------------------

IF ( r(1) == zero ) THEN
  dmrst(1)                   = zero
ELSE
  dmrst(1)                   = frpith * r(1)**3 * rho(2)
END IF

rstmss(1)                    = dmrst(1)

DO j = jr_min,jr_max
  dmrst(j)                   = frpith * ( r(j) * ( r(j) + r(j-1) ) + r(j-1) * r(j-1) ) * ( r(j) - r(j-1) ) * rho(j)
  rstmss(j)                  = rstmss(j-1) + dmrst(j)
END DO

!-----------------------------------------------------------------------
!  Set timestep
!-----------------------------------------------------------------------

dtnph                        = dtime

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute change in psi0 due to neutrino energy advection
!-----------------------------------------------------------------------

CALL nu_energy_advct_x( jr_min, jr_max, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Compute change in psi1 due to neutrino energy advection
!-----------------------------------------------------------------------

CALL psi1_cal( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& rhoa, t, ye, r, rstmss, u, psi0, psi1, nx, nez, nnu, 1 )

!-----------------------------------------------------------------------
!  Compute minimum neutrino energy advection time step
!-----------------------------------------------------------------------

CALL e_advct_time_step( jr_min, jr_max )

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Updated neutrino distribution function
!-----------------------------------------------------------------------

psi0p(imin:imax,:,:,ij_ray,ik_ray)   = psi0(jr_min:jr_max,:,:)
psi1p(imin:imax+1,:,:,ij_ray,ik_ray) = psi1(imin:imax+1,:,:)

!-----------------------------------------------------------------------
!  Time step constraints
!-----------------------------------------------------------------------

dtp (40+1:40+nnu,ij_ray,ik_ray) = dt_nu_e_advct (1:nnu)
jdtp(40+1:40+nnu,ij_ray,ik_ray) = jdt_nu_e_advct(1:nnu)

RETURN
END SUBROUTINE nu_energy_advct_inout_x

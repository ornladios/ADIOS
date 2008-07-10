SUBROUTINE nu_energy_advct_inout_z( imin, imax, kmin, kmax, nx, nz, ki_ray, &
& kj_ray, ij_ray_dim, k_ray_dim, i_radial, j_radial, nez, nnu, nprintp,     &
& rhoip, rhop, rp, rcp, zip, zp, tp, yep, wp, psi0p, psi1p, dtime, rhobar )
!-----------------------------------------------------------------------
!
!    File:         nu_energy_advct_inout_z
!    Module:       nu_energy_advct_inout_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To receive variables from angular_ray_module, perform the neutrino energy,
!       advection step, and return to azimuthal_ray_module the updated variables.
!
!    Input arguments:
!
!  imin              : minimum x-array index
!  imax              : maximum x-array index
!  kmin              : minimum z-array index
!  kmax              : maximum z-array index
!  nx                : x-array extent
!  nz                : z-array extent
!  ki_ray            : x (radial) index of a specific z (azimuthal) ray
!  kj_ray            : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim        : the number of y-zones on a processor before swapping with y
!  k_ray_dim         : the number of radial zones on a processor after swapping with z
!  i_radial          : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  j_radial          : the shifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nez               : neutrino energy array extent
!  nnu               : neutrino flavor array extent
!  nprintp           : unit number to print diagnostics
!  rhoip             : density before hydro advance (cm^{-3})
!  rhop              : density after hydro advance (cm^{-3})
!  rp                : radial zone edge (cm)
!  rcp               : radial zone center (cm)
!  zip               : z (azimuthal) zone position before hydro advance
!  zp                : z (azimuthal) zone position after hydro advance
!  tp                : temperature (K)
!  yep               : electron fraction
!  wp                : z-velocity of zone (cm s^{-1})
!  psi0p             : initial zero angular moments of the neutrino occupation number
!  psi1p             : first angular moments of the neutrino occupation number
!  dtime             : time step
!  rhobar            : mean density at a given radius
!
!    Output arguments:
!
!  psi0p             : updated zero angular moments of the neutrino occupation number
!  psi1p             : first angular moments of the neutrino occupation number
!  dtp               : minimum e-advection time step restrictions (s)
!  jdtp              : zones causing minimum time step restrictions dtp
!
!    Subprograms called:
!  nu_energy_advct_z : performs the neutrino energy advection due to z-hydro
!  e_advct_time_step : computes the neutrino energy advection time step restrictions
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, nu_energy_grid_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith

USE e_advct_module, ONLY : ivc_z, rho, rhoa, z, za, r, rjmh, u, psi0, psi1, &
& jdt_nu_e_advct, dt_nu_e_advct, dpsi0nmax, dtnph, nprint, rhomin_z_eadvect
USE nu_energy_grid_module, ONLY : nnugpmx
USE parallel_module, ONLY : myid
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin             ! minimum x-array index
INTEGER, INTENT(in)              :: imax             ! maximum x-array index
INTEGER, INTENT(in)              :: kmin             ! minimum z-array index
INTEGER, INTENT(in)              :: kmax             ! maximum z-array index
INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: nz               ! z-array extent

INTEGER, INTENT(in)              :: ki_ray           ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray           ! z (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! the number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim        ! the number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)              :: i_radial         ! the unshifted radial zone corresponding to ki_ray, kj_ray
INTEGER, INTENT(in)              :: j_radial         ! the shifted radial zone corresponding to ki_ray, kj_ray

INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent

INTEGER, INTENT(in)              :: nprintp          ! unit number to print diagnostics

REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: rhoip    ! density before Lagrangian hydro step (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: rhop     ! density after Lagrangian hydro step (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                          :: rp       ! radial zone edge (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)                            :: rcp      ! radial zone center (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nz+1)                          :: zip      ! z (azimuthal)
REAL(KIND=double), INTENT(in), DIMENSION(nz+1,ij_ray_dim,k_ray_dim+1)   :: zp       ! z (azimuthal)
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: tp       ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: yep      ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)       :: wp       ! zone-centered z-velocity velocity (cm s^{-1})
REAL(KIND=double), INTENT(in)                                           :: dtime    ! time step
REAL(KIND=double), INTENT(in), DIMENSION(nx)                            :: rhobar   ! mean density at a given radius

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

!INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim, k_ray_dim)                :: jdtp     ! zone causing dtp

!REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim, k_ray_dim)      :: dtp      ! minimum allowed time step

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nz,nez,nnu,ij_ray_dim,k_ray_dim) :: psi0p    ! zero moment of the neutrino occupation probability
REAL(KIND=double), INTENT(inout), DIMENSION(nz,nez,nnu,ij_ray_dim,k_ray_dim) :: psi1p    ! first moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(nz) :: t                ! temperature (K)
REAL(KIND=double), DIMENSION(nz) :: ye               ! electron fraction

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if rhobar < rhomin_z_eadvect
!-----------------------------------------------------------------------

IF ( rhobar( i_radial ) < rhomin_z_eadvect ) RETURN

!-----------------------------------------------------------------------
!  Edit unit number
!-----------------------------------------------------------------------

nprint                     = nprintp

!-----------------------------------------------------------------------
!
!           \\\\\ SET UP FOR NEUTRINO ENERGY ADVECTION /////
!
!  Load variables received from radial_ray_module into e_advct_module
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize time step parameters
!-----------------------------------------------------------------------

dpsi0nmax                  = zero
jdt_nu_e_advct             = 0
dt_nu_e_advct              = 1.d+20

!-----------------------------------------------------------------------
!  Return if nnugpmx = 0  or  ivc_z = 0
!-----------------------------------------------------------------------

IF ( nnugpmx == 0  .or.  ivc_z == 0 ) THEN
!  DO n = 1,nnu
!    dtp(44+n,j_ray)      = dt_nu_e_advct(n)
!    jdtp(44+n,j_ray)     = jdt_nu_e_advct(n)
!  END DO
  RETURN
END IF ! nnugpmx == 0  .or.  ivc_z == 0

!-----------------------------------------------------------------------
!  Transfer psi0 to e_advct_module
!-----------------------------------------------------------------------

psi0(kmin:kmax,:,:)        = psi0p(kmin:kmax,:,:,kj_ray,ki_ray)
psi1(kmin:kmax,:,:)        = psi1p(kmin:kmax,:,:,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Transfer state variables to e_advct_module
!-----------------------------------------------------------------------

rho (kmin:kmax)            = rhoip(kmin:kmax,kj_ray,ki_ray)
rhoa(kmin:kmax)            = rhop (kmin:kmax,kj_ray,ki_ray)
t   (kmin:kmax)            = tp   (kmin:kmax,kj_ray,ki_ray)
ye  (kmin:kmax)            = yep  (kmin:kmax,kj_ray,ki_ray)
u   (kmin:kmax)            = wp   (kmin:kmax,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Transfer x-coordinates to e_advct_module
!-----------------------------------------------------------------------

r   (imin:imax+1)          = rp (imin:imax+1)
rjmh(imin:imax)            = rcp(imin:imax)

!-----------------------------------------------------------------------
!  Transfer z-coordinates to e_advct_module
!-----------------------------------------------------------------------

z (kmin:kmax+1)            = zip(kmin:kmax+1)
za(kmin:kmax+1)            = zp (kmin:kmax+1,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Set timestep
!-----------------------------------------------------------------------

dtnph                      = dtime

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute change in psi0 due to neutrino energy advection
!-----------------------------------------------------------------------

CALL nu_energy_advct_z( kmin, kmax, ki_ray, kj_ray, i_radial )

!-----------------------------------------------------------------------
!  Compute minimum neutrino energy advection time step
!-----------------------------------------------------------------------

CALL e_advct_time_step( kmin, kmax )

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Updated neutrino distribution function
!-----------------------------------------------------------------------

psi0p(kmin:kmax,:,:,kj_ray,ki_ray) = psi0(kmin:kmax,:,:)

!-----------------------------------------------------------------------
!  Time step constraints
!-----------------------------------------------------------------------

!DO n = 1,nnu
!  dtp(44+n,kj_ray,ki_ray)        = dt_nu_e_advct(n)
!  jdtp(44+n,kj_ray,ki_ray)       = jdt_nu_e_advct(n)
!END DO

RETURN
END SUBROUTINE nu_energy_advct_inout_z

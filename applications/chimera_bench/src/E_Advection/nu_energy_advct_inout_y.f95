SUBROUTINE nu_energy_advct_inout_y( imin, imax, jmin, jmax, nx, ny, ji_ray, &
& jk_ray, j_ray_dim, ik_ray_dim, i_radial, j_radial, nez, nnu, nprintp, &
& rhoip, rhop, rp, rcp, yip, yp, tp, yep, vp, psi0p, psi1p, dtime, rhobar )
!-----------------------------------------------------------------------
!
!    File:         nu_energy_advct_inout_y
!    Module:       nu_energy_advct_inout_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To receive variables from angular_ray_module, perform the neutrino energy,
!       advection step, and return to angular_ray_module the updated variables.
!
!    Input arguments:
!
!  imin              : minimum x-array index
!  imax              : maximum x-array index
!  jmin              : minimum y-array index
!  jmax              : maximum y-array index
!  nx                : x-array extent
!  ny                : y-array extent
!  ji_ray            : x (radial) index of a specific y (angular) ray
!  jk_ray            : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim         : the number of radial zones on a processor after swapping with y
!  ik_ray_dim        : the number of z-zones on a processor before swapping with z
!  i_radial          : the unshifted radial zone (angular ray) corresponding to ji_ray, jk_ray
!  j_radial          : the shifted radial zone (angular ray) corresponding to ji_ray, jk_ray
!  nez               : neutrino energy array extent
!  nnu               : neutrino flavor array extent
!  nprintp           : unit number to print diagnostics
!  rhoip             : density before hydro advance (cm^{-3})
!  rhop              : density after hydro advance (cm^{-3})
!  rp                : radial zone edge (cm)
!  rcp               : radial zone center (cm)
!  yip               : angular zone position before hydro advance
!  yp                : angular zone position after hydro advance
!  tp                : temperature (K)
!  yep               : electron fraction
!  vp                : y-velocity of zone (cm s^{-1})
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
!  nu_energy_advct_y : performs the neutrino energy advection due to y-hydro
!  e_advct_time_step : computes the neutrino energy advection time step restrictions
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, nu_energy_grid_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith

USE e_advct_module, ONLY : ivc_y, rho, rhoa, y, ya, r, rjmh, u, psi0, psi1, &
& jdt_nu_e_advct, dt_nu_e_advct, dpsi0nmax, dtnph, nprint, rhomin_y_eadvect
USE nu_energy_grid_module, ONLY : nnugpmx
USE parallel_module, ONLY : myid
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin             ! minimum x-array index
INTEGER, INTENT(in)              :: imax             ! maximum x-array index
INTEGER, INTENT(in)              :: jmin             ! minimum y-array index
INTEGER, INTENT(in)              :: jmax             ! maximum y-array index
INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: ny               ! y-array extent

INTEGER, INTENT(in)              :: ji_ray           ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray           ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: j_ray_dim        ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial         ! the unshifted radial zone corresponding to ji_ray, jk_ray
INTEGER, INTENT(in)              :: j_radial         ! the shifted radial zone corresponding to ji_ray, jk_ray

INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent

INTEGER, INTENT(in)              :: nprintp          ! unit number to print diagnostics

REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: rhoip    ! density before Lagrangian hydro step (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: rhop     ! density after Lagrangian hydro step (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                          :: rp       ! radial zone edge (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)                            :: rcp      ! radial zone center (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny+1)                          :: yip      ! angular zone position before hydro advance
REAL(KIND=double), INTENT(in), DIMENSION(ny+1,j_ray_dim,ik_ray_dim+1)   :: yp       ! angular zone position after hydro advance
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: tp       ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: yep      ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)       :: vp       ! zone-centered y-velocity velocity (cm s^{-1})
REAL(KIND=double), INTENT(in)                                           :: dtime    ! time step
REAL(KIND=double), INTENT(in), DIMENSION(nx)                            :: rhobar   ! mean density at a given radius

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

!INTEGER, INTENT(out), DIMENSION(50,j_ray_dim, ik_ray_dim)                :: jdtp     ! zone causing dtp

!REAL(KIND=double), INTENT(out), DIMENSION(50,j_ray_dim, ik_ray_dim)      :: dtp      ! minimum allowed time step

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim) :: psi0p    ! zero moment of the neutrino occupation probability
REAL(KIND=double), INTENT(inout), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim) :: psi1p    ! first moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(ny) :: t                ! temperature (K)
REAL(KIND=double), DIMENSION(ny) :: ye               ! electron fraction

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if rhobar < rhomin_y_eadvect
!-----------------------------------------------------------------------

IF ( rhobar( i_radial ) < rhomin_y_eadvect ) RETURN

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
!  Return if nnugpmx = 0  or  ivc_y = 0
!-----------------------------------------------------------------------

IF ( nnugpmx == 0  .or.  ivc_y == 0 ) THEN
!  DO n = 1,nnu
!    dtp(44+n,j_ray)      = dt_nu_e_advct(n)
!    jdtp(44+n,j_ray)     = jdt_nu_e_advct(n)
!  END DO
  RETURN
END IF ! nnugpmx == 0  .or.  ivc_y == 0

!-----------------------------------------------------------------------
!  Transfer psi0 to e_advct_module
!-----------------------------------------------------------------------

psi0(jmin:jmax,:,:)        = psi0p(jmin:jmax,:,:,ji_ray,jk_ray)
psi1(jmin:jmax,:,:)        = psi1p(jmin:jmax,:,:,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Transfer state variables to e_advct_module
!-----------------------------------------------------------------------

rho (jmin:jmax)            = rhoip(jmin:jmax,ji_ray,jk_ray)
rhoa(jmin:jmax)            = rhop (jmin:jmax,ji_ray,jk_ray)
t   (jmin:jmax)            = tp   (jmin:jmax,ji_ray,jk_ray)
ye  (jmin:jmax)            = yep  (jmin:jmax,ji_ray,jk_ray)
u   (jmin:jmax)            = vp   (jmin:jmax,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Transfer x-coordinates to e_advct_module
!-----------------------------------------------------------------------

r   (imin:imax+1)          = rp (imin:imax+1)
rjmh(imin:imax)            = rcp(imin:imax)

!-----------------------------------------------------------------------
!  Transfer y-coordinates to e_advct_module
!-----------------------------------------------------------------------

y (jmin:jmax+1)            = yip(jmin:jmax+1)
ya(jmin:jmax+1)            = yp (jmin:jmax+1,ji_ray,jk_ray)

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

CALL nu_energy_advct_y( jmin, jmax, ji_ray, jk_ray, i_radial )

!-----------------------------------------------------------------------
!  Compute minimum neutrino energy advection time step
!-----------------------------------------------------------------------

CALL e_advct_time_step( jmin, jmax )

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Updated neutrino distribution function
!-----------------------------------------------------------------------

psi0p(jmin:jmax,:,:,ji_ray,jk_ray) = psi0(jmin:jmax,:,:)

!-----------------------------------------------------------------------
!  Time step constraints
!-----------------------------------------------------------------------

!DO n = 1,nnu
!  dtp(44+n,ji_ray,jk_ray)        = dt_nu_e_advct(n)
!  jdtp(44+n,ji_ray,jk_ray)       = jdt_nu_e_advct(n)
!END DO

RETURN
END SUBROUTINE nu_energy_advct_inout_y

SUBROUTINE nu_energy_advct_inout_y( imin, imax, jmin, jmax, nx, ny, j_ray, &
& j_ray_dim, i_radial, j_radial, nez, nnu, nprintp, rhoip, rhop, rp, rcp, &
& yip, yp, tp, yep, vp, psi0p, psi1p, dtime, rhobar, nu_strp )
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
!  j_ray             : index denoting a specific angular ray
!  j_ray_dim         : number of angular rays assigned to a processor
!  i_radial          : the unshifted radial zone (angular ray) corresponding to j_ray
!  j_radial          : the shifted radial zone (angular ray) corresponding to j_ray
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
!  nu_strp           : neutrino stress (dynes g^{-1})
!
!    Output arguments:
!
!  psi0p             : updated zero angular moments of the neutrino occupation number
!  psi1p             : first angular moments of the neutrino occupation number
!  nu_strp           : neutrino stress (dynes g^{-1})
!  dtp               : minimum e-advection time step restrictions (s)
!  jdtp              : zones causing minimum time step restrictions dtp
!
!    Subprograms called:
!  nu_energy_advct   : performs the neutrino energy advection
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
USE nu_dist_module, ONLY : stress_y
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
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

INTEGER, INTENT(in)              :: j_ray            ! index denoting a specific angular ray
INTEGER, INTENT(in)              :: j_ray_dim        ! number of angular rays assigned to a processor
INTEGER, INTENT(in)              :: i_radial         ! the unshifted radial zone corresponding to j_ray
INTEGER, INTENT(in)              :: j_radial         ! the shifted radial zone corresponding to j_ray

INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent

INTEGER, INTENT(in)              :: nprintp          ! unit number to print diagnostics

REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim)       :: rhoip    ! density before Lagrangian hydro step (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim)       :: rhop     ! density after Lagrangian hydro step (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)               :: rp       ! radial zone edge (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)                 :: rcp      ! radial zone center (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny+1)               :: yip      ! angular zone position before hydro advance
REAL(KIND=double), INTENT(in), DIMENSION(ny+1,j_ray_dim+1)   :: yp       ! angular zone position after hydro advance
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim)       :: tp       ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim)       :: yep      ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim)       :: vp       ! zone-centered y-velocity velocity (cm s^{-1})
REAL(KIND=double), INTENT(in)                                :: dtime    ! time step
REAL(KIND=double), INTENT(in), DIMENSION(nx)                 :: rhobar   ! mean density at a given radius

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

!INTEGER, INTENT(out), DIMENSION(50,j_ray_dim)                :: jdtp     ! zone causing dtp

!REAL(KIND=double), INTENT(out), DIMENSION(50,j_ray_dim)      :: dtp      ! minimum allowed time step

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(ny,nez,nnu,j_ray_dim) :: psi0p    ! zero moment of the neutrino occupation probability
REAL(KIND=double), INTENT(inout), DIMENSION(ny,nez,nnu,j_ray_dim) :: psi1p    ! first moment of the neutrino occupation probability
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim)         :: nu_strp  ! neutrino stress (dynes g^{-1})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i                ! do index
INTEGER                          :: k                ! neutrino energy index
INTEGER                          :: n                ! neutrino flavor index

REAL(KIND=double), DIMENSION(ny) :: t                ! temperature (K)
REAL(KIND=double), DIMENSION(ny) :: ye               ! electron fraction
REAL(KIND=double), DIMENSION(ny) :: nu_strs_y        ! y-component of neutrino stress

 1001 FORMAT (' Error in closing e_advct_keys.d in subroutine nu_energy_advct_in_y')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Return if rhobar < rhomin_y_eadvect
!-----------------------------------------------------------------------

IF ( rhobar( j_ray + myid * j_ray_dim ) < rhomin_y_eadvect ) RETURN

!-----------------------------------------------------------------------
!        Edit unit number
!-----------------------------------------------------------------------

nprint                   = nprintp

!-----------------------------------------------------------------------
!
!           \\\\\ SET UP FOR NEUTRINO ENERGY ADVECTION /////
!
!        Load variables received from radial_ray_module into
!         e_advct_module
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Initialize time step parameters
!-----------------------------------------------------------------------

print *,' nu_energy_advct_in_y 1, myid=',myid
DO n = 1,nnu
  dpsi0nmax(n)           = zero
  jdt_nu_e_advct(n)      = 0
  dt_nu_e_advct(n)       = 1.d+20
END DO

!-----------------------------------------------------------------------
!        Return if nnugpmx = 0  or  ivc_y = 0
!-----------------------------------------------------------------------

IF ( nnugpmx == 0  .or.  ivc_y == 0 ) THEN
!  DO n = 1,nnu
!    dtp(44+n,j_ray)      = dt_nu_e_advct(n)
!    jdtp(44+n,j_ray)     = jdt_nu_e_advct(n)
!  END DO
  RETURN
END IF

!-----------------------------------------------------------------------
!        Transfer psi0 to e_advct_module
!-----------------------------------------------------------------------

print *,' nu_energy_advct_in_y 2, myid=',myid
DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO i = jmin,jmax
      psi0(i,k,n)        = psi0p(i,k,n,j_ray)
      psi1(i,k,n)        = psi1p(i,k,n,j_ray)
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!        Transfer state variables to e_advct_module
!-----------------------------------------------------------------------

print *,' nu_energy_advct_in_y 3, myid=',myid
DO i = jmin,jmax
  rho(i)                 = rhoip(i,j_ray)
  rhoa(i)                = rhop(i,j_ray)
  t(i)                   = tp(i,j_ray)
  ye(i)                  = yep(i,j_ray)
  u(i)                   = vp(i,j_ray)
END DO

!-----------------------------------------------------------------------
!        Transfer x-coordinates to e_advct_module
!-----------------------------------------------------------------------

print *,' nu_energy_advct_in_y 4, myid=',myid
DO i = imin,imax+1
  r(i)                   = rp(i)
END DO  

DO i = imin,imax
  rjmh(i)                = rcp(i)
END DO  

!-----------------------------------------------------------------------
!        Transfer y-coordinates to e_advct_module
!-----------------------------------------------------------------------

print *,' nu_energy_advct_in_y 5, myid=',myid
DO i = jmin,jmax+1
  y(i)                   = yip(i)
  ya(i)                  = yp(i,j_ray)
END DO

!-----------------------------------------------------------------------
!        Set timestep
!-----------------------------------------------------------------------

dtnph                    = dtime

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO ENERGY ADVECTION STEP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Compute change in psi0 due to neutrino energy advection
!-----------------------------------------------------------------------

print *,' nu_energy_advct_in_y 6, myid=',myid
CALL nu_energy_advct_y( jmin, jmax, j_ray, j_ray_dim )

!-----------------------------------------------------------------------
!        Compute minimum neutrino energy advection time step
!-----------------------------------------------------------------------

print *,' nu_energy_advct_in_y 8, myid=',myid
CALL e_advct_time_step( jmin, jmax )

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Updated neutrino distribution function
!-----------------------------------------------------------------------

print *,' nu_energy_advct_in_y 9, myid=',myid
DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    DO i = jmin,jmax
      psi0p(i,k,n,j_ray) = psi0(i,k,n)
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
!        Time step constraints
!-----------------------------------------------------------------------

!DO n = 1,nnu
!  dtp(44+n,j_ray)        = dt_nu_e_advct(n)
!  jdtp(44+n,j_ray)       = jdt_nu_e_advct(n)
!END DO

RETURN
END SUBROUTINE nu_energy_advct_inout_y

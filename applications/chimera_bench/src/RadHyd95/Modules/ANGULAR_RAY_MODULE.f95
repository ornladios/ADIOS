MODULE angular_ray_module
!-----------------------------------------------------------------------
!
!    File:         angular_ray_module
!    Module:       angular_ray_module
!    Type:         Module
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/01/04
!
!    Purpose:
!      Contains the angular ray variables on a given processor.
!
!    Include files:
!      kind_module
!
!-----------------------------------------------------------------------

USE kind_module
SAVE

!-----------------------------------------------------------------------
!
!           \\\\\ TIME STEPS AND TIME STEP CONSTRAINTS /////
!
!-----------------------------------------------------------------------
!  jdt_y : angular zone setting time step by criterion i, radial zone
!   ji_ray, jk_ray
!  dt_y  : angular ray time step set by criterion i, radial zone ji_ray.
!   jk_ray
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                 :: jdt_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: dt_y
REAL(KIND=double)                                      :: dt_y_state

!-----------------------------------------------------------------------
!
!         \\\\\ PHYSICAL ANGULAR COORDINATES DIFFERENCES /////
!
!-----------------------------------------------------------------------
!  dyphys_c : physical coordinate difference of y_ei(i_1) - uei(i)
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: dyphys_c

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------
!  rho_y     : density: zone average
!  t_y       : temperature: zone average
!  ye_y      : electron fraction: zone average
!  ei_y      : internal energy: zone average
!  p_y       : pressure: zone average
!  gc_y      : 1st adiabatic index
!  ge_y      : 1 + p/(ei*rho)
!  u_y       : zone average velocity x direction
!  v_y       : zone average velocity y direction
!  w_y       : zone average velocity z direction
!  v_e       : zone edged average velocity y-direction
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: rho_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: t_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ye_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ei_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: p_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: gc_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ge_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: u_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: v_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: w_y

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: v_e

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - INITIAL VALUES /////
!
!-----------------------------------------------------------------------
!  rho_yi : density: zone average
!  t_yi   : temperature: zone average
!  ye_yi  : electron fraction: zone average
!  ei_yi  : internal energy: zone average
!  u_yi   : zone average velocity x direction
!  v_yi   : zone average velocity y direction
!  w_yi   : zone average velocity z direction
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: rho_yi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: t_yi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ye_yi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ei_yi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: u_yi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: v_yi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: w_yi

!-----------------------------------------------------------------------
!
!                     \\\\\ GR VARIABLES /////
!
!-----------------------------------------------------------------------
!  agr_y : lapse function
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: agr_y

!-----------------------------------------------------------------------
!
!               \\\\\ SHOCK STABILIZING VARIABLES /////
!
!-----------------------------------------------------------------------
!  flat_y   : variables aligned along angular rays indicating the
!   presence of an angular shock
!
!  flat_x_y : variables aligned along angular rays indicating the
!   presence of a radial shock
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: flat_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: flat_x_y

!-----------------------------------------------------------------------
!
!              \\\\\ GRAVITATION ACCELERATIONS /////
!
!-----------------------------------------------------------------------
!  grav_y_cy  : zone-centered y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_y_ey  : zone-edged y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: grav_y_cy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: grav_y_ey

!-----------------------------------------------------------------------
!
!          \\\\\ COMPOSITION VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------
!  xn_y(j,nc,i)    : mass fraction of the ith nucleus.
!  uburn(j,i)      : cumulative energy generated in zone j by nuclear reactions (ergs/gm).
!  be_nuc_rep(j,i) : binding energy of the representative heavy nucleus (MeV).
!  a_nuc_rep(j,i)  : mass number of the representative heavy nucleus.
!  z_nuc_rep(j,i)  : charge number of the representative heavy nucleus.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)     :: xn_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: uburn_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: be_nuc_rep_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: a_nuc_rep_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: z_nuc_rep_y

!-----------------------------------------------------------------------
!
!                      \\\\\ NSE FLAG /////
!
!-----------------------------------------------------------------------
!  nse_y(j,i) : a nuclear statistical equilibrium flag for angular zone j
!
!     nse_y(j,i) = 0 : material not in nuclear statistical equilibrium;
!      nuclear reaction network must be turned on to evolve the matter
!      composition.
!     nse_y(j,i) = 1 : material in nuclear statistical equilibrium;
!      nuclear reaction network turned off.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                 :: nse_y

!-----------------------------------------------------------------------
!
!                  \\\\\ EOS AND OPACITY GRID /////
!
!-----------------------------------------------------------------------
!  idty_y(j,i) : the density regime (i.e., 1, 2, or 3) of radial zone j
!   as given by the above 
!  inequalities.
!     regime 1:             rho < rhoes(1)
!     regime 2:        rhoes(1) < rho < rhoes(2)
!     regime 3:             rhoes(2) < rho
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                 :: idty_y

!-----------------------------------------------------------------------
!
!                   \\\\\ NEUTRINO ARRAYS /////
!
!-----------------------------------------------------------------------
!  psi0(j,k,n,i,ji_ray,jk_ray) : the zero moment of the occupation
!   distribution for neutrinos at the midpoint of angular zone j, of
!   energy zone k, and of type n.
!
!  psi1(j,k,n,i,ji_ray,jk_ray) : the first moment of the occupation
!   distribution for neutrinos at the outer boundary angular zone j, of
!   energy zone k, and of type n.
!
!  nu_str_cy(j,ji_ray,jk_ray)   : y-component of zone-centered neutrino
!   stress [dynes g^{-1}]
!
!  nu_str_ey(j,ji_ray,jk_ray)   : y-component of zone-edged neutrino
!   stress [dynes g^{-1}]
!
!  e_nu_y(j,ji_ray,jk_ray)     : neutrino energy density [ergs cm^{-3}].
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: psi0_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: psi1_y
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: nu_str_cy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: nu_str_ey
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: e_nu_y

END MODULE angular_ray_module

MODULE azimuthal_ray_module
!-----------------------------------------------------------------------
!
!    File:         azimuthal_ray_module
!    Module:       azimuthal_ray_module
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
!  jdt_z : z (azimuthal) zone setting time step by criterion i, radial
!   zone (ij_ray,ik_ray)
!  dt_z  : z (azimuthal) ray time step set by criterion i, radial zone
!   (ij_ray,ik_ray)
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                 :: jdt_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: dt_z
REAL(KIND=double)                                      :: dt_z_state

!-----------------------------------------------------------------------
!
!              \\\\\ PHYSICAL ANGULAR COORDINATES /////
!
!-----------------------------------------------------------------------
!  dzphys_c : zphys_e(i+1) - zphys_e(i)
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: dzphys_c

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------
!  rho_z     : density: zone average
!  t_z       : temperature: zone average
!  ye_z      : electron fraction: zone average
!  ei_z      : internal energy: zone average
!  p_z       : pressure: zone average
!  gc_z      : 1st adiabatic index
!  ge_z      : 1 + p/(ei*rho)
!  u_z       : zone average velocity x direction
!  v_z       : zone average velocity y direction
!  w_z       : zone average velocity z direction
!  w_e       : zone edged average velocity z-direction
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: rho_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: t_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ye_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ei_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: p_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: gc_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ge_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: u_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: v_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: w_z

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: w_e

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - INITIAL VALUES /////
!
!-----------------------------------------------------------------------
!  rho_zi : density: zone average
!  t_zi   : temperature: zone average
!  ye_zi  : electron fraction: zone average
!  ei_zi  : internal energy: zone average
!  u_zi   : zone average velocity x direction
!  v_zi   : zone average velocity y direction
!  w_zi   : zone average velocity z direction
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: rho_zi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: t_zi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ye_zi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: ei_zi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: u_zi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: v_zi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: w_zi

!-----------------------------------------------------------------------
!
!                     \\\\\ GR VARIABLES /////
!
!-----------------------------------------------------------------------
!  agr_z : lapse function
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: agr_z

!-----------------------------------------------------------------------
!
!               \\\\\ SHOCK STABILIZING VARIABLES /////
!
!-----------------------------------------------------------------------
!  flat_z   : variables aligned along azimuthal rays indicating the
!   presence of an azimuthal shock
!
!  flat_x_z : variables aligned along azimuthal rays iindicating the
!   presence of a radial shock
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: flat_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: flat_x_z

!-----------------------------------------------------------------------
!
!              \\\\\ GRAVITATION ACCELERATIONS /////
!
!-----------------------------------------------------------------------
!  grav_z_cz  : zone-centered z-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_z_ez zone-edged z-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: grav_z_cz
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: grav_z_ez

!-----------------------------------------------------------------------
!
!          \\\\\ COMPOSITION VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------
!  xn_z(j,nc,i)    : mass fraction of the ith nucleus.
!  uburn(j,i)      : cumulative energy generated in zone j by nuclear
!   reactions (ergs/gm).
!  be_nuc_rep(j,i) : binding energy of the representative heavy nucleus
!   [MeV].
!  a_nuc_rep(j,i)  : mass number of the representative heavy nucleus.
!  z_nuc_rep(j,i)  : charge number of the representative heavy nucleus.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)     :: xn_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: uburn_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: be_nuc_rep_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: a_nuc_rep_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: z_nuc_rep_z

!-----------------------------------------------------------------------
!
!                      \\\\\ NSE FLAG /////
!
!-----------------------------------------------------------------------
!  nse_z(j,i) : a nuclear statistical equilibrium flag for angular zone j
!
!     nse_z(j,i) = 0 : material not in nuclear statistical equilibrium;
!      nuclear reaction network must be turned on to evolve the matter
!      composition.
!     nse_z(j,i) = 1 : material in nuclear statistical equilibrium;
!      nuclear reaction network turned off.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                 :: nse_z

!-----------------------------------------------------------------------
!
!                  \\\\\ EOS AND OPACITY GRID /////
!
!-----------------------------------------------------------------------
!  idty_z(j,i) : the density regime (i.e., 1, 2, or 3) of radial zone j as
!   given by the above 
!  inequalities.
!     regime 1:             rho < rhoes(1)
!     regime 2:        rhoes(1) < rho < rhoes(2)
!     regime 3:             rhoes(2) < rho
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                 :: idty_z

!-----------------------------------------------------------------------
!
!                   \\\\\ NEUTRINO ARRAYS /////
!
!-----------------------------------------------------------------------
!  psi0(j,k,n,i,ki_ray,kj_ray) : the zero moment of the occupation
!   distribution for neutrinos at the midpoint of angular zone j, of
!   energy zone k, and of type n.
!
!  psi1(j,k,n,i,ki_ray,kj_ray) : the first moment of the occupation
!   distribution for neutrinos at the outer boundary angular zone j,
!   of energy zone k, and of type n.
!
!  nu_str_cy(j,ki_ray,kj_ray)   : y-component of zone-centered neutrino
!   stress [dynes g^{-1}]
!
!  nu_str_ey(j,ki_ray,kj_ray)   : y-component of zone-edged neutrino
!   stress [dynes g^{-1}]
!
!  e_nu_z(j,ki_ray,kj_ray)     : neutrino energy density [ergs cm^{-3}].
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: psi0_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: psi1_z
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: nu_str_cz
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: nu_str_ez
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)       :: e_nu_z

END MODULE azimuthal_ray_module

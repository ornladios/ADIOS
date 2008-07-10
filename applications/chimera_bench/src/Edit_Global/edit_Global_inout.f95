SUBROUTINE edit_Global_inout( imin, imax, nx, nez, nnu, jmin, jmax, &
& ij_ray_dim, ny, kmin, kmax, ik_ray_dim, nz, x_e_in, x_c_in, dx_in, &
& uMD_in, vMD_in, wMD_in, rhoMD_in, tMD_in, yeMD_in, psi0_in, psi1_in, &
& unu_in, dunu_in, unue_in, dunue_in, dtnph_in, time_in, t_bounce_in, &
& ncycle_in, nnc, xn_in, nse_in, be_nuc_rep_in, a_nuc_rep_in, z_nuc_rep_in, &
& d_omega )
!-----------------------------------------------------------------------
!
!    File:         edit_Global_inout
!    Module:       edit_Global_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To test criteria and load variables for the global edits.
!
!    Subprograms called:
!
!  edit_MD_exec : executes the MD edits
!
!    Input arguments:
!
!  imin          : minimum x-array index for the edit
!  imax          : maximum x-array index for the edit
!  nx            : x-array extent
!  nez           : neutrino energy array extent
!  nnu           : neutrino flavor array extent
!  jmin          : minimum y-array index for the edit
!  jmax          : maximum y-array index for the edit
!  ij_ray_dim    : number of y-zones on a processor before swapping with y
!  ny            : y_array extent
!  kmin          : minimum z-array index for the edit
!  kmax          : maximum z-array index for the edit
!  ik_ray_dim    : number of z-zones on a processor before swapping with z
!  nz            : z_array extent
!  x_e_in        : radial edge of zone (cm)
!  x_c_in        : radial midpoint of zone (cm)
!  dx_in         : radial zone width (cm)
!  uMD_in        : radial velocity of zone (cm s^{-1})
!  vMD_in        : angular velocity of zone (cm s^{-1})
!  wMD_in        : angular velocity of zone (cm s^{-1})
!  rhoMD_in      : density of zone (g cm^{-3})
!  tMD_in        : temperature of zone (K)
!  yeMD_in       : electron fraction of zone
!  psi0_in       : zero moment of the neutrino distribution
!  psi1_in       : first moment of the neutrino distribution
!  unu_in        : radial zone-centered neutrino energy
!  dunu_in       : radial zone-centered neutrino energy zone width
!  unue_in       : radial zone-edged neutrino energy
!  dunue_in      : radial zone-edged neutrino energy zone width
!  dtnph_in      : hydro time step
!  time_in       : elapsed time
!  t_bounce_in   : time of core bounce
!  ncycle_in     : cycle number
!  nnc           : composition array extent
!  xn_in         : composition mass fractions
!  nse_in        : NSE flag
!  be_nuc_rep_in : binding energy of auxiliary nucleus (MeV)
!  a_nuc_rep_in  : mass number of auxiliary nucleus
!  z_nuc_rep_in  : charge number of auxiliary nucleus
!  d_omega       : solid angles subtended by radial rays
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  cycle_module, edit_module, eos_snc_x_module, nucbrn_module,
!  nu_dist_module, shock_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc
USE numerical_module, ONLY : zero, third, half, frpi, epsilon
USE physcnst_module, ONLY : msolar, ergfoe, kmev

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nned_global, inted_global, ied_global_n, nprint, &
& nlog, dt_global_ed1, dt_global_ed2, ied_global_t
USE eos_snc_x_module, ONLY : aesv, nuc_number
USE nucbrn_module, ONLY: a_name, dudt_nuc, a_nuc
USE nu_dist_module, ONLY : dudt_nu
USE shock_module, ONLY: j_shk_radial_all_p
USE t_cntrl_module, ONLY: dtnph, time, t_bounce

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                             :: imin             ! minimum x-array index for the edit
INTEGER, INTENT(in)                             :: imax             ! maximum x-array index for the edit
INTEGER, INTENT(in)                             :: nx               ! x-array extent
INTEGER, INTENT(in)                             :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)                             :: nnu              ! neutrino flavor array extent

INTEGER, INTENT(in)                             :: jmin             ! minimum y-array index for the edit
INTEGER, INTENT(in)                             :: jmax             ! number of radial rays assigned to a processor
INTEGER, INTENT(in)                             :: ij_ray_dim       ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)                             :: ny               ! y_array extent

INTEGER, INTENT(in)                             :: kmin             ! minimum y-array index for the edit
INTEGER, INTENT(in)                             :: kmax             ! maximum z-array index for the edit
INTEGER, INTENT(in)                             :: ik_ray_dim       ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)                             :: nz               ! z_array extent

INTEGER, INTENT(in)                             :: ncycle_in        ! ncycle
INTEGER, INTENT(in)                             :: nnc              ! composition array extent
INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: nse_in  ! NSE flag

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)  :: x_e_in           ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: x_c_in           ! radial edge of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: dx_in            ! radial zone width (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: uMD_in    ! radial velocity of zone (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: vMD_in    ! y (angular) velocity of zone (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: wMD_in    ! z (azimuthal) velocity of zone (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: rhoMD_in  ! density of zone (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: tMD_in    ! temperature of zone (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: yeMD_in   ! entropy of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0_in  ! zero moment of the neutrino distribution
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1_in  ! first moment of the neutrino distribution
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: unu_in       ! radial zone-centered neutrino energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: dunu_in      ! radial zone-centered neutrino energy zone width
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: unue_in      ! radial zone-edged neutrino energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: dunue_in     ! radial zone-edged neutrino energy zone width
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim) :: xn_in        ! composition mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: be_nuc_rep_in    ! binding energy of auxiliary nucleus (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: a_nuc_rep_in     ! mass number of auxiliary nucleus
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: z_nuc_rep_in     ! charge number of auxiliary nucleus
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz) :: d_omega          ! solid angles subtended by radial rays

REAL(KIND=double), INTENT(in)                   :: dtnph_in         ! time step used at current cycle
REAL(KIND=double), INTENT(in)                   :: time_in          ! elapsed time
REAL(KIND=double), INTENT(in)                   :: t_bounce_in      ! time of core bounce

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=2), DIMENSION(nx)                :: c_shock          ! shock location

LOGICAL                                         :: first = .true.
LOGICAL                                         :: l_global_ned     ! edit_global flag from cycle criteria
LOGICAL                                         :: l_global_ted     ! edit_global flag from time criteria

INTEGER                                         :: i                ! radial zone index
INTEGER                                         :: j                ! y (angular) zone index
INTEGER                                         :: k                ! z (azimuthal) zone index
INTEGER                                         :: l                ! composition index
INTEGER                                         :: n                ! neutrino energy index
INTEGER                                         :: itime            ! used to determine when to generate a MD edit
INTEGER, SAVE                                   :: itimeprev        ! used to determine when to generate a MD edit

INTEGER, PARAMETER                              :: nleg = 10        ! highest Legendre polynomial to be used 
INTEGER, PARAMETER                              :: ncd = 300        ! composition dimension
INTEGER                                         :: ii               ! composition index
INTEGER                                         :: n_nucp1          ! nuc_number + 1
INTEGER                                         :: n_hvy            ! number of heavy nuclei (not counting representative heavy nucleus)
INTEGER, SAVE                                   :: i_Ni             ! nickel index
INTEGER, SAVE                                   :: i_He             ! helium index
INTEGER, SAVE                                   :: i_neut           ! neutron index
INTEGER, SAVE                                   :: i_prot           ! proton index
INTEGER, SAVE, DIMENSION(ncd)                   :: i_hvy            ! heavy nucleus index

INTEGER, DIMENSION(nx,ny,nz)                    :: nseMD            ! NSE flag
INTEGER                                         :: nse_min          ! minimum value of NSE-nonNSE boundary
INTEGER                                         :: nse_max          ! maximum value of NSE-nonNSE boundary

INTEGER                                         :: j_shock_min      ! minimum radial index of shock
INTEGER                                         :: j_shock_max      ! maximum radial index of shock

REAL(KIND=double)                               :: t_tb             ! time from bounce (s)
REAL(KIND=double)                               :: tmult            ! used to determine when to generate a MD edit

REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: sMD_in    ! entropy of zone
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: pMD_in    ! pressure of zone (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: gamMD_in  ! adiabatic exponent of zone
REAL(KIND=double), DIMENSION(nx)                :: mach_2           ! mach number squared
REAL(KIND=double), DIMENSION(nx)                :: mach_min         ! minimum mach number/radius
REAL(KIND=double), DIMENSION(nx)                :: mach_max         ! maximum mach number/radius
REAL(KIND=double), DIMENSION(nx)                :: mach_mean        ! mass averaged mach number/radius
REAL(KIND=double), DIMENSION(nx)                :: sigma_mach       ! RMS mach number/radius
REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: lum_in           ! neutrino luminosity (foes)
REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: e_rms_stat_in    ! SQRT( SUM psi0 * w5dw/SUM w3ww ) (MeV)
REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: e_rms_trns_in    ! SQRT( SUM psi1 * w5dw/SUM w3ww ) (MeV)

REAL(KIND=double), DIMENSION(nx,ny,nz)          :: volume_zone      ! volume/zone (cm^{3})
REAL(KIND=double), DIMENSION(nx)                :: volume_a_ray     ! volume/(angular ray) (cm^{3})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: mass_zone        ! mass/zone (g)
REAL(KIND=double), DIMENSION(nx)                :: mass_a_ray       ! mass/angular ray (g)
REAL(KIND=double), DIMENSION(nx)                :: rho_mean         ! mass averaged density/radius (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)                :: rho_min          ! minimum density/(radial zone) (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)                :: rho_max          ! maximum density/(radial zone) (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)                :: sigma_rho        ! RMS density/(radial zone) (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)                :: t_mean           ! mass averaged temperature/radius (MeV)
REAL(KIND=double), DIMENSION(nx)                :: t_min            ! minimum temperature/(radial zone) (MeV)
REAL(KIND=double), DIMENSION(nx)                :: t_max            ! maximum temperature/(radial zone) (MeV)
REAL(KIND=double), DIMENSION(nx)                :: sigma_t          ! RMS temperature/(radial zone) (MeV)
REAL(KIND=double), DIMENSION(nx)                :: u_mean           ! mass averaged radial velocity/radius (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: u_min            ! minimum radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: u_max            ! maximum radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: sigma_u          ! RMS radial velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: v_mean           ! mass averaged y-velocity/radius (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: v_min            ! minimum y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: v_max            ! maximum y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: sigma_v          ! RMS y-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: w_mean           ! mass averaged z-velocity/radius (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: w_min            ! minimum z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: w_max            ! maximum z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: sigma_w          ! RMS z-velocity/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: s_mean           ! mass averaged entropy/(radial zone)
REAL(KIND=double), DIMENSION(nx)                :: s_min            ! minimum entropy/(radial zone)
REAL(KIND=double), DIMENSION(nx)                :: s_max            ! maximum entropy/(radial zone)
REAL(KIND=double), DIMENSION(nx)                :: sigma_s          ! RMS entropy/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: ye_mean          ! mass averaged electron fraction/(radial zone)
REAL(KIND=double), DIMENSION(nx)                :: ye_min           ! minimum electron fraction/(radial zone)
REAL(KIND=double), DIMENSION(nx)                :: ye_max           ! maximum electron fraction/(radial zone)
REAL(KIND=double), DIMENSION(nx)                :: sigma_ye         ! RMS electron fraction/(radial zone) (cm s^{-1})
REAL(KIND=double), DIMENSION(nx)                :: dudt_nuc_mean    ! mass averaged energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx)                :: dudt_nuc_min     ! minimum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx)                :: dudt_nuc_max     ! maximum energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx)                :: dudt_nuc_zone    ! energy generation rate by nuclear reactions (B s^{-1} radial-zone^{-1})
REAL(KIND=double), DIMENSION(nx)                :: dudt_nu_mean     ! mass averaged neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx)                :: dudt_nu_min      ! minimum neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx)                :: dudt_nu_max      ! maximum neutrino energy deposition rate (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx)                :: dudt_nu_zone     ! neutrino energy deposition rate  (B s^{-1} radial-zone^{-1})
REAL(KIND=double), DIMENSION(nx,nnu)            :: lum_mean         ! mass averaged n-neutrino luminosity (foes s^{1})
REAL(KIND=double), DIMENSION(nx,nnu)            :: lum_min          ! minimum n-neutrino luminosity (foes s^{1})
REAL(KIND=double), DIMENSION(nx,nnu)            :: lum_max          ! maximum  n-neutrino luminosity (foes s^{1})
REAL(KIND=double), DIMENSION(nx,nnu)            :: sigma_lum        ! RMS neutrino n-neutrino luminosity (foes s^{1})
REAL(KIND=double), DIMENSION(nx,nnu)            :: e_rms_mean       ! mass averaged neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), DIMENSION(nx,nnu)            :: e_rms_min        ! minimum neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), DIMENSION(nx,nnu)            :: e_rms_max        ! maximum neutrino n-neutrino rms energy (MeV)
REAL(KIND=double), DIMENSION(nx,nnu)            :: sigma_e_rms      ! RMS neutrino n-neutrino rms energy (MeV))

REAL(KIND=double), DIMENSION(nx,ny,nz)          :: uMD              ! radial velocity of zone (cm s^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: vMD              ! y (angular) velocity of zone (cm s^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: wMD              ! z (azimuthal) velocity of zone (cm s^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: rhoMD            ! density (g cm^{-3})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: tMD              ! temperature (NeV)
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: sMD              ! entropy of zone
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: yeMD             ! electron fraction of zone
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: pMD              ! pressure of zone (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: gamMD            ! adiabatic exponent of zone
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: machMD           ! mach number
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: e_bindMD         ! binding energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: e_no_bindMD      ! total minus binding energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: e_bind_fnlMD     ! final binding energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: a_nuc_repMD      ! mass number of the representative heavy nucleus
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: dudt_nucMD       ! energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: dudt_nuMD        ! energy deposition rate by neutrinos (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx,nnu,ny,nz)      :: lumMD            ! neutrino luminosity (foes)
REAL(KIND=double), DIMENSION(nx,nnu,ny,nz)      :: e_rmsMD          ! SQRT( SUM psi0 * w5dw/SUM w3ww ) (MeV)

REAL(KIND=double), DIMENSION(nx,ny,nz)          :: grav_pot         ! gravitational potential energy  (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)                :: grav_pot_a_ray   ! angular ray gravitational potential (ergs)
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: ke               ! kinetic energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)                :: ke_a_ray         ! angular ray kinetic energy (ergs)
REAL(KIND=double), DIMENSION(nnc)               :: xn_ijk           ! composition of zone ij
REAL(KIND=double)                               :: e_ph             ! photon energy (ergs g^{-3})
REAL(KIND=double)                               :: e_elec           ! electron energy (ergs g^{-3})
REAL(KIND=double)                               :: e_drip           ! drip energy (ergs g^{-1})
REAL(KIND=double)                               :: e_hvy            ! nuclear internal energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: e_bind_in        ! binding energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: e_no_bind_in     ! total minus binding energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: e_bind_fnl_in    ! final binding energy (ergs g^{-1})
REAL(KIND=double)                               :: e_no_bind        ! total minus binding energy (ergs g^{-1})
REAL(KIND=double)                               :: e_total          ! total internal energy (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx)                :: e_bind_a_ray     ! angular ray current binding energy (ergs)
REAL(KIND=double), DIMENSION(nx)                :: e_no_bind_a_ray  ! angular ray internal minus binding energy (ergs)
REAL(KIND=double), DIMENSION(nx)                :: e_bind_fnl_a_ray ! angular ray final binding energy (ergs)

REAL(KIND=double), DIMENSION(nx,nnc,ny,nz)      :: xnMD             ! composition mass fractions
REAL(KIND=double), DIMENSION(nx,nnc)            :: xn_a_ray         ! composition by mass (solar masses)
REAL(KIND=double), DIMENSION(nx,nnc)            :: xn_mean          ! angular average of the composition mass fractions
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: a_mean_zone      ! angular average of the mean nuclear mass number for each zone
REAL(KIND=double), DIMENSION(nx)                :: a_mean           ! angular average of the mean nuclear mass number along an spherical shell
REAL(KIND=double), DIMENSION(nx)                :: p_mean           ! volume averaged pressure (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nleg)              :: p_SASI           ! power in the convection or SASI as a function of the Legendre mode

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Need MPI version of edit_Global_inout since n_proc=',i4,' > 1')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Stop if n_proc > 1.
!-----------------------------------------------------------------------

IF ( n_proc > 1 ) THEN
  WRITE (nprint,101) n_proc
  WRITE (nlog,101) n_proc
  STOP
END IF ! n_proc > 1

!-----------------------------------------------------------------------
!
!              \\\\\ LOAD VARIABLES INTO MODULES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer integer scalar parameters
!-----------------------------------------------------------------------

ncycle                       = ncycle_in

!-----------------------------------------------------------------------
!  Transfer real scalar parameters
!-----------------------------------------------------------------------

dtnph                        = dtnph_in
time                         = time_in
t_bounce                     = t_bounce_in
t_tb                         = zero

!-----------------------------------------------------------------------
!
!               \\\\\ DETERMINE WHETHER TO EDIT /////
!               /////  -CYCLE NUMBER CRITERIA-  \\\\\
!
!-----------------------------------------------------------------------

l_global_ned                 = .false.

IF ( ied_global_n /= 0 ) THEN
  nned_global                = nned_global + 1
  IF ( nned_global >= inted_global ) THEN
    nned_global              = 0
    l_global_ned             = .true.
  END IF ! nned_global >= inted_global
END IF ! ied_global_n /= 0

!-----------------------------------------------------------------------
!
!               \\\\\ DETERMINE WHETHER TO EDIT /////
!               /////      -TIME CRITERIA-      \\\\\
!
!-----------------------------------------------------------------------

l_global_ted                 = .false.

IF ( ied_global_t /= 0 ) THEN

  IF ( t_bounce <= zero ) THEN

!-----------------------------------------------------------------------
!  Initialize time criterion
!-----------------------------------------------------------------------

    IF ( first ) THEN
      tmult                  = 1.d+3/dt_global_ed1
      itimeprev              = int( time * tmult )
      first                  = .false.
    END IF

!-----------------------------------------------------------------------
!  Generate a global edit every multiple of dt_global_ed1 ms up to
!   bounce
!-----------------------------------------------------------------------

    tmult                    = 1.d+3/dt_global_ed1
    itime                    = int( time * tmult )
    IF ( itime /= itimeprev ) THEN
      l_global_ted           = .true.
      itimeprev              = itime
    END IF
    
  ELSE ! t_bounce > zero

    t_tb                     = time - t_bounce

!-----------------------------------------------------------------------
!  Generate a global edit at bounce
!-----------------------------------------------------------------------

    IF ( t_tb - dtnph <= zero  .and.  t_tb > zero ) l_global_ted = .true.

!-----------------------------------------------------------------------
!  Generate a global edit 1 ms after bounce
!-----------------------------------------------------------------------

    IF ( t_tb - dtnph <= 1.d-3  .and.  t_tb > 1.d-3 ) l_global_ted = .true.

!-----------------------------------------------------------------------
!  Generate a global edit 10 ms after bounce
!-----------------------------------------------------------------------

    IF ( t_tb - dtnph <= 1.d-2  .and.  t_tb > 1.d-2 ) l_global_ted = .true.

!-----------------------------------------------------------------------
!  Initialize time criterion
!-----------------------------------------------------------------------

    IF ( first ) THEN
      tmult                  = 1.d+3/dt_global_ed2
      itimeprev              = int( t_tb * tmult )
      first                  = .false.
    END IF

!-----------------------------------------------------------------------
!  Generate a global edit every multiple of dt_global_ed2 ms after
!   bounce
!-----------------------------------------------------------------------

    tmult                    = 1.d+3/dt_global_ed2
    itime                    = int( t_tb * tmult )
    IF ( itime /= itimeprev ) THEN
      l_global_ted           = .true.
      itimeprev              = itime
    END IF

  END IF ! t_bounce <= zero

END IF ! ied_global_t /= 0

!-----------------------------------------------------------------------
!  Return if all edit_Global criteria false
!-----------------------------------------------------------------------

IF (       ( .not. l_global_ned  ) &
&   .and.  ( .not. l_global_ted  ) ) RETURN

!-----------------------------------------------------------------------
!
!                \\\\\                           /////
!                      COMPUTE VARIABLES TO EDIT
!                /////                           \\\\\
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Zone volumes and masses
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    volume_zone(imin:imax,j,k)                                            &
&                            = ( dx_in(imin:imax) * ( x_e_in(imin:imax) * x_e_in(imin+1:imax+1) &
&                            + dx_in(imin:imax) * dx_in(imin:imax) * third ) ) * d_omega(j,k)
  END DO ! j
END DO ! k

DO i = imin,imax
  volume_a_ray(i)            = SUM( volume_zone(i,:,:) )
END DO

rhoMD(imin:imax,jmin:jmax,kmin:kmax)                                      &
&                            = rhoMD_in(imin:imax,jmin:jmax,kmin:kmax)
mass_zone(imin:imax,jmin:jmax,kmin:kmax)                                  &
&                            = volume_zone(imin:imax,jmin:jmax,kmin:kmax) &
&                            * rhoMD(imin:imax,jmin:jmax,kmin:kmax)
mass_zone(imax+1,jmin:jmax,kmin:kmax)                                     &
&                            = mass_zone(imax,jmin:jmax,kmin:kmax)

DO i = imin,imax
  mass_a_ray(i)              = SUM( mass_zone(i,:,:) )
END DO
mass_a_ray(imax+1)           = mass_a_ray(imax)

!-----------------------------------------------------------------------
!  Composition indexing
!-----------------------------------------------------------------------

n_hvy                        = 0
i_He                         = 0
i_neut                       = 0
i_prot                       = 0
i_hvy                        = 0
ii                           = 0
n_nucp1                      = nuc_number + 1

DO i = 1,nuc_number
  IF (      a_name(i) == '  n  ' ) THEN
    i_neut                   = i
  ELSE IF ( a_name(i) == '  p  ' ) THEN
    i_prot                   = i
  ELSE IF ( a_name(i) == '  4He' ) THEN
    i_He                     = i
  ELSE
    ii                       = ii + 1
    n_hvy                    = n_hvy + 1
    i_hvy(ii)                = i
  END IF
END DO ! i

DO i = 1,nuc_number
  IF ( a_name(i) == ' 56Ni' ) THEN
    i_Ni                     = i
  END IF
END DO ! i

!-----------------------------------------------------------------------
!
!             \\\\\ MEAN, MINIMUM, MAXIMUM, AND RMS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Minimum and maximum radial shock index
!-----------------------------------------------------------------------

j_shock_min                  = MAX( MINVAL(j_shk_radial_all_p(:,:)), 1 )
j_shock_max                  = MAX( MAXVAL(j_shk_radial_all_p(:,:)), 1 )
c_shock                      = '  '
IF ( j_shock_min == 1 ) THEN
  c_shock(j_shock_min)       = '  '
ELSE IF ( j_shock_min == j_shock_max ) THEN
  c_shock(j_shock_min)       = '**'
ELSE
  c_shock(j_shock_min)       = '* '
  c_shock(j_shock_max)       = ' *'
END IF

!-----------------------------------------------------------------------
!  Minimum, maximum, volume averaged, and RMS densities
!-----------------------------------------------------------------------

DO i = imin,imax
  rho_min(i)                 = MINVAL( rhoMD(i,:,:) )
  rho_max(i)                 = MAXVAL( rhoMD(i,:,:) )
  rho_mean(i)                = SUM( rhoMD(i,jmin:jmax,kmin:kmax)                              &
&                            * volume_zone(i,jmin:jmax,kmin:kmax) )/volume_a_ray(i)
  sigma_rho(i)               = DSQRT( SUM( ( rhoMD(i,jmin:jmax,kmin:kmax) -  rho_mean(i) )**2 &
&                            * volume_zone(i,jmin:jmax,kmin:kmax) )/volume_a_ray(i) + epsilon )
END DO ! i

!-----------------------------------------------------------------------
!  Minimum, maximum, mass averaged, and RMS temperatures
!-----------------------------------------------------------------------

tMD(imin:imax,jmin:jmax,kmin:kmax)                                                      &
&                            = kmev * tMD_in(imin:imax,jmin:jmax,kmin:kmax)
DO i = imin,imax
  t_min(i)                   = MINVAL( tMD(i,:,:) )
  t_max(i)                   = MAXVAL( tMD(i,:,:) )
  t_mean(i)                  = SUM( tMD(i,jmin:jmax,kmin:kmax)                            &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  sigma_t(i)                 = DSQRT( SUM( ( tMD(i,jmin:jmax,kmin:kmax) -  t_mean(i) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
END DO ! i

!-----------------------------------------------------------------------
!  Minimum, maximum, mass averaged, and RMS x-velocities
!-----------------------------------------------------------------------

uMD(imin:imax,jmin:jmax,kmin:kmax)                                                      &
&                            = uMD_in(imin:imax,jmin:jmax,kmin:kmax)
DO i = imin,imax
  u_min(i)                   = MINVAL( uMD(i,:,:) )
  u_max(i)                   = MAXVAL( uMD(i,:,:) )
  u_mean(i)                  = SUM( uMD(i,jmin:jmax,kmin:kmax)                            &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  sigma_u(i)                 = DSQRT( SUM( ( uMD(i,jmin:jmax,kmin:kmax) -  u_mean(i) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
END DO ! i

!-----------------------------------------------------------------------
!  Minimum, maximum, mass averaged, and RMS y-velocities
!-----------------------------------------------------------------------

vMD(imin:imax,jmin:jmax,kmin:kmax)                                                      &
&                            = vMD_in(imin:imax,jmin:jmax,kmin:kmax)
DO i = imin,imax
  v_min(i)                   = MINVAL( vMD(i,:,:) )
  v_max(i)                   = MAXVAL( vMD(i,:,:) )
  v_mean(i)                  = SUM( vMD(i,jmin:jmax,kmin:kmax)                            &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  sigma_v(i)                 = DSQRT( SUM( ( vMD(i,jmin:jmax,kmin:kmax) -  v_mean(i) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
END DO ! i

!-----------------------------------------------------------------------
!  Minimum, maximum, mass averaged, and RMS z-velocities
!-----------------------------------------------------------------------

wMD(imin:imax,jmin:jmax,kmin:kmax)                                                      &
&                            = wMD_in(imin:imax,jmin:jmax,kmin:kmax)
DO i = imin,imax
  w_min(i)                   = MINVAL( wMD(i,:,:) )
  w_max(i)                   = MAXVAL( wMD(i,:,:) )
  w_mean(i)                  = SUM( wMD(i,jmin:jmax,kmin:kmax)                            &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  sigma_w(i)                 = DSQRT( SUM( ( wMD(i,jmin:jmax,kmin:kmax) -  w_mean(i) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
END DO ! i

!-----------------------------------------------------------------------
!  Minimum, maximum, mass averaged, and RMS entropies
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    sMD_in(imin:imax,j,k)    = aesv(imin+1:imax+1,3,j,k)
  END DO ! j
END DO ! j

sMD(imin:imax,jmin:jmax,kmin:kmax)                                                      &
&                            = sMD_in(imin:imax,jmin:jmax,kmin:kmax)
DO i = imin,imax
  s_min(i)                   = MINVAL( sMD(i,:,:) )
  s_max(i)                   = MAXVAL( sMD(i,:,:) )
  s_mean(i)                  = SUM( sMD(i,jmin:jmax,kmin:kmax)                            &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  sigma_s(i)                 = DSQRT( SUM( ( sMD(i,jmin:jmax,kmin:kmax) -  s_mean(i) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
END DO ! i

!-----------------------------------------------------------------------
!  Minimum, maximum, mass averaged, and RMS electron fractions
!-----------------------------------------------------------------------

yeMD(imin:imax,jmin:jmax,kmin:kmax)                                                       &
&                            = yeMD_in(imin:imax,jmin:jmax,kmin:kmax)
DO i = imin,imax
  ye_min(i)                  = MINVAL( yeMD(i,:,:) )
  ye_max(i)                  = MAXVAL( yeMD(i,:,:) )
  ye_mean(i)                 = SUM( yeMD(i,jmin:jmax,kmin:kmax)                             &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  sigma_ye(i)                = DSQRT( SUM( ( yeMD(i,jmin:jmax,kmin:kmax) -  ye_mean(i) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
END DO ! i

!-----------------------------------------------------------------------
!  Mass averaged mach number
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    pMD_in  (imin:imax,j,k)  = aesv(imin+1:imax+1,1,j,k)
    gamMD_in(imin:imax,j,k)  = aesv(imin+1:imax+1,12,j,k)
  END DO ! j
END DO ! j

pMD  (imin:imax,jmin:jmax,kmin:kmax)                                    &
&                            = pMD_in  (imin:imax,jmin:jmax,kmin:kmax)
gamMD(imin:imax,jmin:jmax,kmin:kmax)                                    &
&                            = gamMD_in(imin:imax,jmin:jmax,kmin:kmax)
DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    mach_2(imin:imax)        = ( uMD(imin:imax,j,k) * uMD(imin:imax,j,k)  &
&                            + vMD(imin:imax,j,k) * vMD(imin:imax,j,k) + epsilon ) &
&                            * rhoMD(imin:imax,j,k)/( gamMD(imin:imax,j,k) * pMD(imin:imax,j,k) )
    machMD(imin:imax,j,k)    = DSQRT( mach_2(imin:imax) )
  END DO ! j
END DO ! k

DO i = imin,imax
  mach_min(i)                = MINVAL( machMD(i,:,:) )
  mach_max(i)                = MAXVAL( machMD(i,:,:) )
  mach_mean(i)               = SUM( machMD(i,jmin:jmax,kmin:kmax)                               &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  sigma_mach(i)              = DSQRT( SUM( ( machMD(i,jmin:jmax,kmin:kmax) -  mach_mean(i) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
END DO ! i

!-----------------------------------------------------------------------
!
!          \\\\\ NEUTRINO LUMINOSITIES AND RMS ENERGIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute the neutrino luminosities
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    CALL luminosity_MD( imin, imax, j, k, ij_ray_dim, ik_ray_dim, nx, nez, &
&    nnu, x_e_in, psi1_in, unue_in, dunue_in, lum_in )
  END DO ! j
END DO ! k

!-----------------------------------------------------------------------
!  Minimum, maximum, mass averaged, and RMS neutrino luminosities
!-----------------------------------------------------------------------

lumMD(imin:imax,:,jmin:jmax,kmin:kmax)                                  &
&                            = lum_in(imin:imax,:,jmin:jmax,kmin:kmax)
DO n = 1,nnu
    DO i = imin,imax+1
    lum_min(i,n)             = MINVAL( lumMD(i,n,:,:) )
    lum_max(i,n)             = MAXVAL( lumMD(i,n,:,:) )
    lum_mean(i,n)            = SUM( lumMD(i,n,jmin:jmax,kmin:kmax)                              &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
    sigma_lum(i,n)           = DSQRT( SUM( ( lumMD(i,n,jmin:jmax,kmin:kmax) -  lum_mean(i,n) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
  END DO ! i
END DO ! n

!-----------------------------------------------------------------------
!  Compute neutrino rms energies
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    CALL e_rms_MD( imin, imax, j, k, ij_ray_dim, ik_ray_dim, nx, nez, nnu, &
&    psi0_in, psi1_in, unu_in, dunu_in, e_rms_stat_in, e_rms_trns_in )
  END DO ! j
END DO ! k

!-----------------------------------------------------------------------
!  Minimum, maximum, mass averaged, and RMS neutrino RMS energies
!-----------------------------------------------------------------------

e_rmsMD(imin:imax,:,jmin:jmax,kmin:kmax)                                &
&                            = e_rms_stat_in(imin:imax,:,jmin:jmax,kmin:kmax)
DO n = 1,nnu
  DO i = imin,imax
    e_rms_min(i,n)           = MINVAL( e_rmsMD(i,n,:,:) )
    e_rms_max(i,n)           = MAXVAL( e_rmsMD(i,n,:,:) )
    e_rms_mean(i,n)          = SUM( e_rmsMD(i,n,jmin:jmax,kmin:kmax)                                &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
    sigma_e_rms(i,n)         = DSQRT( SUM( ( e_rmsMD(i,n,jmin:jmax,kmin:kmax) -  e_rms_mean(i,n) )**2 &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i) + epsilon )
  END DO ! i
END DO ! n

!-----------------------------------------------------------------------
!
!                       \\\\\ ENERGETICS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Minimum, maximum, and mass averaged energy generation rates by
!   nuclear reactions
!-----------------------------------------------------------------------

dudt_nucMD(imin:imax,jmin:jmax,kmin:kmax)                                                     &
&                          = dudt_nuc(imin+1:imax+1,jmin:jmax,kmin:kmax)
DO i = imin,imax
  dudt_nuc_min(i)          = MINVAL( dudt_nucMD(i,:,:) )
  dudt_nuc_max(i)          = MAXVAL( dudt_nucMD(i,:,:) )
  dudt_nuc_mean(i)         = SUM( dudt_nucMD(i,jmin:jmax,kmin:kmax)                           &
&                          * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  dudt_nuc_zone(i)         = dudt_nuc_mean(i) * mass_a_ray(i) * ergfoe
END DO ! i

!-----------------------------------------------------------------------
!  Minimum, maximum, and mass averaged neutrino energy deposition rates
!-----------------------------------------------------------------------

dudt_nuMD(imin:imax,jmin:jmax,kmin:kmax)                                                      &
&                          = dudt_nu(imin+1:imax+1,jmin:jmax,kmin:kmax)
DO i = imin,imax
  dudt_nu_min(i)           = MINVAL( dudt_nuMD(i,:,:) )
  dudt_nu_max(i)           = MAXVAL( dudt_nuMD(i,:,:) )
  dudt_nu_mean(i)          = SUM( dudt_nuMD(i,jmin:jmax,kmin:kmax)                            &
&                          * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  dudt_nu_zone(i)          = dudt_nu_mean(i) * mass_a_ray(i) * ergfoe
END DO ! i

!-----------------------------------------------------------------------
!  Gravitational potential
!-----------------------------------------------------------------------

CALL grav_potential( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, x_e_in, &
& mass_a_ray, grav_pot )

DO i = imin,imax
  grav_pot_a_ray(i)        = SUM( grav_pot(i,jmin:jmax,kmin:kmax)                             &
&                          * mass_zone(i,jmin:jmax,kmin:kmax) )
END DO ! i

!-----------------------------------------------------------------------
!  Kinetic energy
!-----------------------------------------------------------------------

ke(imin:imax,jmin:jmax,kmin:kmax)                                         &
&                          = half * ( uMD(imin:imax,jmin:jmax,kmin:kmax)  &
&                          * uMD(imin:imax,jmin:jmax,kmin:kmax)           &
&                          + vMD(imin:imax,jmin:jmax,kmin:kmax)           &
&                          * vMD(imin:imax,jmin:jmax,kmin:kmax)           &
&                          + wMD(imin:imax,jmin:jmax,kmin:kmax)           &
&                          * wMD(imin:imax,jmin:jmax,kmin:kmax) )
DO i = imin,imax
  ke_a_ray(i)              = SUM( ke(i,jmin:jmax,kmin:kmax) * mass_zone(i,jmin:jmax,kmin:kmax) )
END DO ! i

!-----------------------------------------------------------------------
!  Internal energy and nuclear binding energy
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    DO i = imin,imax
      xn_ijk(1:nnc)        = xn_in(i,1:nnc,j,k)
      CALL eos_nnse_e( i+1, rhoMD_in(i,j,k), tMD_in(i,j,k), yeMD_in(i,j,k), xn_ijk, &
&      nnc, a_nuc_rep_in(i,j,k), z_nuc_rep_in(i,j,k), be_nuc_rep_in(i,j,k), e_ph,   &
&      e_elec, e_drip, e_hvy, e_bind_in(i,j,k), e_no_bind_in(i,j,k), e_total )
      e_no_bind_in(i,j,k)  = e_ph + e_elec + e_drip + e_hvy
      xn_ijk               = zero
      xn_ijk(i_Ni)         = 1.d0 - xn_in(i,n_nucp1,j,k)
      CALL eos_nnse_e( j+1, rhoMD_in(i,j,k), tMD_in(i,j,k), yeMD_in(i,j,k), xn_ijk, &
&      nnc, a_nuc_rep_in(i,j,k), z_nuc_rep_in(i,j,k), be_nuc_rep_in(i,j,k), e_ph,   &
&      e_elec, e_drip, e_hvy, e_bind_fnl_in(i,j,k), e_no_bind, e_total )
    END DO ! i
  END DO ! j
END DO ! k

e_no_bindMD(imin:imax,jmin:jmax,kmin:kmax)                              &
&                            = e_no_bind_in(imin:imax,jmin:jmax,kmin:kmax)
e_bindMD(imin:imax,jmin:jmax,kmin:kmax)                                 &
&                            = e_bind_in(imin:imax,jmin:jmax,kmin:kmax)
e_bind_fnlMD(imin:imax,jmin:jmax,kmin:kmax)                             &
&                            = e_bind_fnl_in(imin:imax,jmin:jmax,kmin:kmax)

!-----------------------------------------------------------------------
!  Internal energy minus the binding energy
!-----------------------------------------------------------------------

DO i = imin,imax
  e_no_bind_a_ray(i)         = SUM( e_no_bindMD(i,jmin:jmax,kmin:kmax)  &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )
END DO ! i

!-----------------------------------------------------------------------
!  Current binding energy
!-----------------------------------------------------------------------

DO i = imin,imax
  e_bind_a_ray(i)            = SUM( e_bindMD(i,jmin:jmax,kmin:kmax)     &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )
END DO ! i

!-----------------------------------------------------------------------
!  Final binding energy
!-----------------------------------------------------------------------

DO i = imin,imax
  e_bind_fnl_a_ray(i)        = SUM( e_bind_fnlMD(i,jmin:jmax,kmin:kmax) &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )
END DO ! 

!-----------------------------------------------------------------------
!
!                       \\\\\ COMPOSITION /////
!
!-----------------------------------------------------------------------

xnMD(imin:imax,1:nnc,jmin:jmax,kmin:kmax)                               &
&                            = xn_in(imin:imax,1:nnc,jmin:jmax,kmin:kmax)
nseMD(imin:imax,jmin:jmax,kmin:kmax)                                    &
&                            = nse_in(imin:imax,jmin:jmax,kmin:kmax) 
a_nuc_repmD(imin:imax,jmin:jmax,kmin:kmax)                              &
&                            = a_nuc_rep_in(imin:imax,jmin:jmax,kmin:kmax)

!-----------------------------------------------------------------------
!  Minimum and maximum NSE-nonNSE boundary
!-----------------------------------------------------------------------

nse_min                      = nx + 1
nse_max                      = nx + 1

DO i = imin,imax
  IF ( MINVAL(nseMD(i,:,:)) == 0 ) THEN
  nse_min                    = i
  EXIT
  END IF ! MINVAL(nseMD(i,:)) == 0
END DO ! i = imin,imax


DO i = imin,imax
  IF ( MAXVAL(nseMD(i,:,:)) == 0 ) THEN
  nse_max                    = i
  EXIT
  END IF ! MAXVAL(nseMD(i,:)) == 0
END DO ! i = imin,imax

!-----------------------------------------------------------------------
!  Composition masses
!-----------------------------------------------------------------------

DO l = 1,nnc
  DO i = imin,imax
    xn_a_ray(i,l)            = SUM( xnMD(i,l,jmin:jmax,kmin:kmax)       &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/msolar
  END DO ! i = imin,imax
END DO ! l = 1,nnc

!-----------------------------------------------------------------------
!  Composition mass fractions
!-----------------------------------------------------------------------

DO l = 1,nnc
  DO i = imin,imax
    xn_mean(i,l)             = xn_a_ray(i,l) * msolar/mass_a_ray(i)
  END DO ! i = imin,imax
END DO ! l = 1,nnc

!-----------------------------------------------------------------------
!  Mean nuclear mass number
!-----------------------------------------------------------------------

  DO i = imin,imax
    DO j = jmin,jmax
      DO k = kmin,kmax
        a_mean_zone(i,j,k)   = SUM( xnMD(i,1:nuc_number,j,k) * a_nuc(1:nuc_number) )
        a_mean_zone(i,j,k)   = a_mean_zone(i,j,k) + a_nuc_repMD(i,j,k) * xnMD(i,n_nucp1,j,k)
      END DO ! k = kmin,kmax
    END DO ! j = jmin,jmax
  END DO ! i = imin,imax

  DO i = imin,imax
    a_mean(i)                = SUM( a_mean_zone(i,jmin:jmax,kmin:kmax)  &
&                            * mass_zone(i,jmin:jmax,kmin:kmax) )/mass_a_ray(i)
  END DO ! i = imin,imax

!-----------------------------------------------------------------------
!
!    \\\\\ CONVECTION AND SASI POWER AS A FUNCTION OF L-MODE  /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Mean pressure
!-----------------------------------------------------------------------

  DO i = imin,imax
    p_mean(i)                = SUM( pMD(i,jmin:jmax,kmin:kmax)          &
&                            * volume_zone(i,jmin:jmax,kmin:kmax) )/volume_a_ray(i)
  END DO ! i = imin,imax

!-----------------------------------------------------------------------
!  Power as a function of l-mode
!-----------------------------------------------------------------------

!  CALL power_SASI( imin, imax, nx, ny, jsphere_min - 1, j_shock_max - 1, &
!&  x_e_in, x_c_in, dx_in, y_in, pMD, p_mean, p_SASI )

!-----------------------------------------------------------------------
!
!               \\\\\ TRANSFER VARIABLES TO EDIT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Execute the MD-edits
!-----------------------------------------------------------------------

CALL edit_Global_exec( imin, imax, nx, nnu, time, t_tb, x_c_in, mass_a_ray, &
& rho_min, rho_max, rho_mean, sigma_rho, t_min, t_max, t_mean, sigma_t, &
& u_min, u_max, u_mean, sigma_u, v_min, v_max, v_mean, sigma_v, w_min, &
& w_max, w_mean, sigma_w, s_min, s_max, s_mean, sigma_s, ye_min, ye_max, &
& ye_mean, sigma_ye, dudt_nuc_min, dudt_nuc_max, dudt_nuc_mean, dudt_nuc_zone, &
& dudt_nu_min, dudt_nu_max, dudt_nu_mean, dudt_nu_zone, mach_min, mach_max, &
& mach_mean, sigma_mach, lum_min, lum_max, lum_mean, sigma_lum, e_rms_min, &
& e_rms_max, e_rms_mean, sigma_e_rms, grav_pot_a_ray, ke_a_ray, e_no_bind_a_ray, &
&  e_bind_a_ray, e_bind_fnl_a_ray, xn_a_ray, xn_mean, a_mean, nnc, nse_min, &
&  nse_max, c_shock, p_SASI, l_global_ned, l_global_ted )

RETURN
END SUBROUTINE edit_Global_inout

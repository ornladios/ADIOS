SUBROUTINE edit_MD_inout( imin, imax, nx, nez, nnu, jmin, jmax, ij_ray_dim, &
& ny, kmin, kmax, ik_ray_dim, nz, x_e_in, x_c_in, y_in, z_in, uMD_in,       &
& vMD_in, wMD_in, rhoMD_in, tMD_in, yeMD_in, rhobar, psi0_in, psi1_in,      &
& e_nu_MD_in, f_nu_MD_in, unu_in, dunu_in, unue_in, dunue_in, dtnph_in,     &
& time_in, t_bounce_in, xn_in, nse_in, nnc, ncycle_in, grav_x_e_in,         &
& grav_x_c_in, grav_y_c_in, grav_z_c_in, gtot_pot_c )
!-----------------------------------------------------------------------
!
!    File:         edit_MD_inout_MPI
!    Module:       edit_MD_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To test criteria and load variables for the MD edits.
!
!    Subprograms called:
!
!  edit_MD_exec : executes the MD edits
!
!    Input arguments:
!
!  imin         : minimum x-array index for the edit
!  imax         : maximum x-array index for the edit
!  nx           : x-array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  jmin         : minimum y-array index for the edit
!  jmax         : maximum y-array index for the edit
!  ij_ray_dim   : number of y-zones on a processor before swapping with y
!  ny           : y-array extent
!  kmin         : minimum z-array index for the edit
!  kmax         : maximum z-array index for the edit
!  ik_ray_dim   : number of z-zones on a processor before swapping with z
!  nz           : z-array extent
!  x_e_in       : radial edge of zone [cm]
!  x_c_in       : radial midpoint of zone [cm]
!  y_in         : y (angular) midpoint of zone
!  z_in         : z (azimuthal) midpoint of zone
!  uMD_in       : x (radial) velocity of zone [cm s^{-1}}
!  vMD_in       : y (angular) velocity of zone [cm s^{-1}}
!  wMD_in       : z (azimuthal) velocity of zone [cm s^{-1}}
!  rhoMD_in     : density of zone [g cm^{-3}]
!  tMD_in       : temperature of zone [K]
!  yeMD_in      : electron fraction of zone
!  rhobar       : mean density at a given radius [g cm^{-3}]
!  psi0_in      : zero moment of the neutrino distribution
!  psi1_in      : first moment of the neutrino distribution
!  e_nu_MD_in   : neutrino energy density [ergs cm^{-3}]
!  f_nu_MD_in   : neutrino energy flux [ergs cm^{-2} s^{-1}]
!  unu_in       : radial zone-centered neutrino energy
!  dunu_in      : radial zone-centered neutrino energy zone width
!  unue_in      : radial zone-edged neutrino energy
!  dunue_in     : radial zone-edged neutrino energy zone width
!  dtnph_in     : hydro time step
!  time_in      : elapsed time
!  t_bounce_in  : time of core bounce
!  xn_in        : composition mass fractions
!  nse_in       : NSE flag
!  nnc          : composition array extent
!  ncycle_in    : cycle number
!  grav_x_e_in  : zone-edged x-component of gravitational acceleration [cm s^{-2} g^{-1}]
!  grav_x_c_in  : zone-centered x-component of gravitational acceleration [cm s^{-2} g^{-1}]
!  grav_y_c_in  : zone-centered y-component of gravitational acceleration [cm s^{-2} g^{-1}]
!  grav_z_c_in  : zone-centered z-component of gravitational acceleration [cm s^{-2} g^{-1}]
!  gtot_pot_c   : unshifted zone-centered gravitational potential energy [ergs]
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, eos_snc_x_module, nucbrn_module,
!  nu_dist_module, nu_energy_grid_module, parallel_module, shock_module,
!  t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc, n_proc_y, n_proc_z
USE numerical_module, ONLY : zero, third, half, one, epsilon, frpi, frpith, &
& ncoef
USE physcnst_module, ONLY : ergmev, rmu

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nedMDu, nedMDv, nedMDw, nedMDs, nedMDd, nedMDe,          &
& nedMDp, nedMDenu, nedMDfnu, nedMDa, nedMDx, nedMDye, nedMDcm, nedMDnu,         &
& nedMDnc, nedMDnl, nedMDne, nedMDgx, nedMDgy, nedMDgz, nedMDBVw, nedMDyl,       &
& intedMDu, intedMDv, intedMDw, intedMDs, intedMDd, intedMDe, intedMDp,          &
& intedMDenu, intedMDfnu, intedMDa, intedMDx, intedMDye, intedMDcm, intedMDnu,   &
& intedMDnc, intedMDnl, intedMDne, intedMDgx, intedMDgy, intedMDgz, intedMDBVw,  &
& intedMDyl, iedMDu, iedMDv, iedMDw, iedMDs, iedMDd, iedMDe, iedMDp, iedMDenu,   &
& iedMDfnu, iedMDa, iedMDgy, iedMDgz, iedMDBVw, iedMDyl, dt_MDedit1, dt_MDedit2, &
& iedMDx, iedMDye, iedMDcm, iedMDnu, iedMDnc, iedMDnl, iedMDne, iedMDgx,         &
& i_editMD, nlog
USE eos_snc_x_module, ONLY : aesv, nuc_number
USE nucbrn_module, ONLY: a_name, dudt_nuc, a_nuc, a_nuc_rep
USE nu_dist_module, ONLY : dunujeadt, dudt_nu
USE nu_energy_grid_module, ONLY : nnugp
USE parallel_module, ONLY : myid, ierr
USE shock_module, ONLY : pq_x, pqy_x, j_shk_radial_p
USE t_cntrl_module, ONLY: dtnph, time, t_bounce

USE mpi

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                             :: imin            ! minimum x-array index for the edit
INTEGER, INTENT(in)                             :: imax            ! maximum x-array index for the edit
INTEGER, INTENT(in)                             :: nx              ! x-array extent
INTEGER, INTENT(in)                             :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)                             :: nnu             ! neutrino flavor array extent

INTEGER, INTENT(in)                             :: jmin            ! minimum y-array index for the edit
INTEGER, INTENT(in)                             :: jmax            ! number of radial rays assigned to a processor
INTEGER, INTENT(in)                             :: ij_ray_dim      ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)                             :: ny              ! y_array extent

INTEGER, INTENT(in)                             :: kmin            ! minimum y-array index for the edit
INTEGER, INTENT(in)                             :: kmax            ! maximum z-array index for the edit
INTEGER, INTENT(in)                             :: ik_ray_dim      ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)                             :: nz              ! z_array extent

INTEGER, INTENT(in)                             :: ncycle_in       ! ncycle
INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: nse_in ! NSE flag
INTEGER, INTENT(in)                             :: nnc             ! composition array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)  :: x_e_in          ! radial midpoint of zone [cm]
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: x_c_in          ! radial edge of zone [cm]
REAL(KIND=double), INTENT(in), DIMENSION(ny)    :: y_in            ! y (angular) midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nz)    :: z_in            ! z (azimuthal) midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: uMD_in   ! x (radial) velocity of zone [cm s^{-1}}
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: vMD_in   ! y (angular) velocity of zone [cm s^{-1}}
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: wMD_in   ! z (azimuthal) velocity of zone [cm s^{-1}}
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: rhoMD_in ! density of zone [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: tMD_in   ! temperature of zone [K]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: yeMD_in  ! entropy of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rhobar          ! mean density at a given radius [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0_in  ! zero moment of the neutrino distribution
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1_in  ! first moment of the neutrino distribution
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: e_nu_MD_in ! neutrino energy density [ergs cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)       :: f_nu_MD_in ! neutrino energy flux (ergs s^{-1} cm^{-2})
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: unu_in       ! radial zone-centered neutrino energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: dunu_in      ! radial zone-centered neutrino energy zone width
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: unue_in      ! radial zone-edged neutrino energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim) :: dunue_in     ! radial zone-edged neutrino energy zone width
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim) :: xn_in        ! composition mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: grav_x_e_in  ! zone-edged x-component of gravitational acceleration [cm s^{-2} g^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: grav_x_c_in  ! zone-centered x-component of gravitational acceleration [cm s^{-2} g^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: grav_y_c_in  ! zone-centered y-component of gravitational acceleration [cm s^{-2} g^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: grav_z_c_in  ! zone-centered z-component of gravitational acceleration [cm s^{-2} g^{-1}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: gtot_pot_c   ! unshifted zone-centered gravitational potential energy [ergs g^{-1}]

REAL(KIND=double), INTENT(in)                   :: dtnph_in        ! time step used at current cycle
REAL(KIND=double), INTENT(in)                   :: time_in         ! elapsed time
REAL(KIND=double), INTENT(in)                   :: t_bounce_in     ! time of core bounce

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=14)                              :: var_name

LOGICAL                                         :: first = .true.
LOGICAL                                         :: first_c = .true.
LOGICAL                                         :: l_editMDu       ! editMDu flag from cycle criteria
LOGICAL                                         :: l_editMDv       ! editMDv flag from cycle criteria
LOGICAL                                         :: l_editMDw       ! editMDw flag from cycle criteria
LOGICAL                                         :: l_editMDs       ! editMDs flag from cycle criteria
LOGICAL                                         :: l_editMDd       ! editMDd flag from cycle criteria
LOGICAL                                         :: l_editMDe       ! editMDe flag from cycle criteria
LOGICAL                                         :: l_editMDp       ! editMDp flag from cycle criteria
LOGICAL                                         :: l_editMDenu     ! editMDenu flag from cycle criteria
LOGICAL                                         :: l_editMDfnu     ! editMDfnu flag from cycle criteria
LOGICAL                                         :: l_editMDa       ! editMDa flag from cycle criteria
LOGICAL                                         :: l_editMDx       ! editMDx flag from cycle criteria
LOGICAL                                         :: l_editMDye      ! editMDye flag from cycle criteria
LOGICAL                                         :: l_editMDcm      ! editMDcm flag from cycle criteria
LOGICAL                                         :: l_editMDnu      ! editMDnu flag from cycle criteria
LOGICAL                                         :: l_editMDnc      ! editMDnc flag from cycle criteria
LOGICAL                                         :: l_editMDnl      ! editMDnl flag from cycle criteria
LOGICAL                                         :: l_editMDne      ! editMDne flag from cycle criteria
LOGICAL                                         :: l_editMDgx      ! editMDgx flag from cycle criteria
LOGICAL                                         :: l_editMDgy      ! editMDgy flag from cycle criteria
LOGICAL                                         :: l_editMDgz      ! editMDgz flag from cycle criteria
LOGICAL                                         :: l_editMDBVw     ! editMDBVw flag from cycle criteria
LOGICAL                                         :: l_editMDyl      ! editMDyl flag from cycle criteria
LOGICAL                                         :: l_MDedit        ! editMD  flag from time criteria

INTEGER                                         :: i               ! radial zone index
INTEGER                                         :: jr              ! shifted radial zone index
INTEGER                                         :: j               ! y (angular) zone index
INTEGER                                         :: k               ! z (azimuthal) zone index
INTEGER                                         :: n               ! neutrino energy index
INTEGER                                         :: m               ! processor index
INTEGER                                         :: mj              ! y-block index
INTEGER                                         :: mk              ! z-block index
INTEGER                                         :: jsk             ! y-array index of gathered array
INTEGER                                         :: ksk             ! z-array index of gathered array

INTEGER                                         :: istat           ! allocation status

INTEGER                                         :: jr_min          ! shifted miminim radial zone index
INTEGER                                         :: jr_max          ! shifted maximum radial zone index
INTEGER                                         :: jr_maxp         ! jr_max + 1
INTEGER                                         :: nwt             ! standard deviation flag for linear fitting

INTEGER                                         :: ic              ! EOS index
INTEGER, DIMENSION(ny,nz)                       :: i_nse           ! radial zone index of NSE-nonNSE boundary
INTEGER                                         :: i_16O           ! 16O abundance index
INTEGER                                         :: i_O1u           ! x_16O > 0.1 index
INTEGER                                         :: i_O1l           ! x_16O < 0.1 index
INTEGER                                         :: i_xOu           ! x_16O > x_O index
INTEGER                                         :: i_xOl           ! x_16O < x_O index
INTEGER                                         :: i_12C           ! 12C abundance index
INTEGER                                         :: i_20Ne          ! 20Ne abundance index
INTEGER                                         :: i_24Mg          ! 24Mg abundance index
INTEGER                                         :: i_28Si          ! 28Si abundance index
INTEGER                                         :: i_32S           ! 28Si abundance index
INTEGER                                         :: i_36Ar          ! 36Ar abundance index
INTEGER                                         :: i_40Ca          ! 40Ca abundance index
INTEGER                                         :: i_44Ti          ! 44Ti abundance index
INTEGER                                         :: i_48Cr          ! 48Cr abundance index
INTEGER                                         :: i_52Fe          ! 52Fe abundance index
INTEGER                                         :: i_56Ni          ! 56Ni abundance index
INTEGER                                         :: i_60Zn          ! 60Zn abundance index
INTEGER                                         :: i_n             ! neutron abundance index
INTEGER                                         :: i_p             ! proton abundance index
INTEGER                                         :: i_4He           ! 4He abundance index
INTEGER                                         :: j_shock         ! shifted radial zone index of shock maximum
INTEGER                                         :: j_shock_mn      ! shifted radial zone index of minimum estimated shock radius
INTEGER                                         :: j_shock_mx      ! shifted radial zone index of maximum estimated shock radius
INTEGER                                         :: i_shock         ! radial zone index of shock maximum
INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim)       :: i_shock_in      ! radial zone index of shock maximum
INTEGER                                         :: i_shockm        ! radial zone index of shock next to maximum
INTEGER                                         :: i_shock_mn      ! radial zone index of minimum estimated shock radius
INTEGER                                         :: i_shock_mx      ! radial zone index of maximum estimated shock radius
INTEGER                                         :: itime           ! used to determine when to generate a MD edit
INTEGER                                         :: itimeprev       ! used to determine when to generate a MD edit
INTEGER, DIMENSION(nx,ny,nz)                    :: nseMD           ! NSE flag
INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim)       :: j_shock_in      ! shifted radial zone index of shock maximum
INTEGER, DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: j_sphere        ! zone index of the neutrinospheres
INTEGER, DIMENSION(nnu,ij_ray_dim,ik_ray_dim)   :: jsphere_mean_in ! shifted radial zone <= mean n-neutrinosphere
INTEGER                                         :: j_gain          ! shifted radial index of gain radius
INTEGER                                         :: i_gain          ! unshifted radial index of gain radius
INTEGER, DIMENSION(nx,ny,nz)                    :: i_comp          ! integer denoting the dominant element
INTEGER, DIMENSION(1)                           :: i_comp_mx       ! index of the dominant element of r_comp

INTEGER                                         :: c_gath_recv     ! gather recv buffer count
INTEGER                                         :: c_gath_send     ! gather send buffer count

INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim)       :: i_gain_in       ! index of the gain radius

INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:)         :: irecv_buf      ! receive buffer for gathering data from processors

REAL(KIND=double), PARAMETER                    :: UTOT0 = 8.9d0   ! change in the zero of energy (MeV)
REAL(KIND=double), PARAMETER                    :: ku = ergmev/rmu ! ( # nucleons/gram )( erg/mev )
REAL(KIND=double), PARAMETER                    :: e_bar0   = ku * UTOT0

REAL(KIND=double), DIMENSION(nx)                :: rho             ! shifted density [cm^{-3}]
REAL(KIND=double), DIMENSION(nx)                :: t               ! shifted temperature [K]
REAL(KIND=double), DIMENSION(nx)                :: ye              ! electron fraction
REAL(KIND=double), DIMENSION(nx)                :: ye_yl           ! electron - lepton fraction
REAL(KIND=double), DIMENSION(nx)                :: g_acc_e         ! adged component of the gravitational acceleration [cm s^{-2} g^{-1}]
REAL(KIND=double), DIMENSION(nx)                :: r               ! enclosed rest mass [cm]
REAL(KIND=double), DIMENSION(nx)                :: rstmss          ! enclosed rest mass [g]

REAL(KIND=double)                               :: t_tb            ! time from bounce [s]
REAL(KIND=double)                               :: tmult           ! used to determine when to generate a MD edit
REAL(KIND=double)                               :: v_csoundMD_2    ! square of the mach number
REAL(KIND=double), PARAMETER                    :: pqmin = 1.d0    ! shock criterion

REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: pMD_in          ! pressure of zone [ergs cm^{-3}]
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: eMD_in          ! internal energy of zone (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: sMD_in          ! entropy of zone
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: aMD_in          ! mean nuclear mass number of zone
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: v_csoundMD_in   ! mach number
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: e_nud_MD_in     ! neutrino energy per unit mass (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: lum_in          ! neutrino luminosity (foes)
REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: e_rms_stat_in   ! SQRT( SUM psi0 * w5dw/SUM w3ww ) (MeV)
REAL(KIND=double), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: e_rms_trns_in   ! SQRT( SUM psi1 * w5dw/SUM w3ww ) (MeV)
REAL(KIND=double)                               :: q_shock_in1     ! strength of shock maximum
REAL(KIND=double)                               :: q_shock_in2     ! strength of shock second to maximum
REAL(KIND=double)                               :: r_shock_in1     ! radius of shock maximum
REAL(KIND=double)                               :: r_shock_in2     ! radius of shock second to maximum
REAL(KIND=double)                               :: r_shk           ! radius of shock maximum
REAL(KIND=double)                               :: m_shk           ! mass enclosed by shock maximum
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim)         :: r_shock_in      ! radius of shock maximum
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim)         :: r_shock_in_mn   ! minimum estimateed shock radius
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim)         :: r_shock_in_mx   ! maximum estimateed shock radius
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)     :: rsphere_mean_in ! mean neutrinosphere radius
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)     :: dsphere_mean_in ! mean neutrinosphere density
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)     :: tsphere_mean_in ! mean neutrinosphere temperature
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)     :: msphere_mean_in ! mean neutrinosphere enclosed mass
REAL(KIND=double), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)     :: esphere_mean_in ! mean neutrinosphere energy
REAL(KIND=double), DIMENSION(nnu+1,ij_ray_dim,ik_ray_dim)   :: r_gain_in       ! gain radius
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim)         :: tau_adv_in      ! advection time scale (s)
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim)         :: tau_heat_nu_in  ! neutrino heating time scale (s)
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim)         :: tau_heat_nuc_in ! nuclear heating time scale (s)
REAL(KIND=double), DIMENSION(ij_ray_dim,ik_ray_dim)         :: n_grow_in       ! number of convective e-foldings in fluid from shock to gain radius

REAL(KIND=double), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: r_sphere        ! neutrinosphere radius
REAL(KIND=double), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: d_sphere        ! neutrinosphere density
REAL(KIND=double), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: t_sphere        ! neutrinosphere temperature
REAL(KIND=double), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: m_sphere        ! neutrinosphere enclosed mass

REAL(KIND=double), DIMENSION(nx,ny,nz)          :: rhoMD           ! density [g cm^{-3}]
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: uMD             ! x (radial) velocity of zone [cm s^{-1}}
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: vMD             ! y (angular) velocity of zone [cm s^{-1}}
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: wMD             ! z (azimuthal) velocity of zone [cm s^{-1}}
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: pMD             ! pressure of zone [ergs cm^{-3}]
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: eMD             ! internal energy of zone (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: sMD             ! entropy of zone
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: aMD             ! mean nuclear mass number of zone
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: yeMD            ! electron fraction of zone
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: v_csoundMD      ! mach number
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: e_nu_MD         ! neutrino energy density (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx+1,ny,nz)        :: f_nu_MD         ! neutrino flux (ergs s^{-1} cm^{-2})
REAL(KIND=double), DIMENSION(nx,nnu,ny,nz)      :: lumMD           ! neutrino luminosity (foes)
REAL(KIND=double), DIMENSION(nx,nnu,ny,nz)      :: e_rmsMD         ! SQRT( SUM psi0 * w5dw/SUM w3ww ) (MeV)
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: grav_x_cMD      ! zone-centered x-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: grav_y_cMD      ! zone-centered y-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: grav_z_cMD      ! zone-centered z-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), DIMENSION(nx,nnc,ny,nz)      :: xnMD            ! composition mass fractions
REAL(KIND=double), DIMENSION(nx,12,ny,nz)       :: aesvMD          ! EOS variables

REAL(KIND=double), DIMENSION(ny,nz)             :: r_nse           ! radius of NSE-nonNSE boundary
REAL(KIND=double), DIMENSION(ny,nz)             :: r_O1            ! radius of x(16O)=0.1 boundary
REAL(KIND=double), DIMENSION(ny,nz)             :: r_xO            ! radius of x(16O)=x_O boundary
REAL(KIND=double), PARAMETER                    :: x_O = 0.3d0
REAL(KIND=double), DIMENSION(ny,nz)             :: r_shock         ! radius of shock maximum
REAL(KIND=double), DIMENSION(ny,nz)             :: r_shock_mn      ! minimum estimateed shock radius
REAL(KIND=double), DIMENSION(ny,nz)             :: r_shock_mx      ! maximum estimateed shock radius
REAL(KIND=double), DIMENSION(nnu,ny,nz)         :: rsphere_mean    ! mean neutrinosphere radius
REAL(KIND=double), DIMENSION(nnu,ny,nz)         :: dsphere_mean    ! mean neutrinosphere density
REAL(KIND=double), DIMENSION(nnu,ny,nz)         :: tsphere_mean    ! mean neutrinosphere temperature
REAL(KIND=double), DIMENSION(nnu,ny,nz)         :: msphere_mean    ! mean neutrinosphere enclosed mass
REAL(KIND=double), DIMENSION(nnu,ny,nz)         :: esphere_mean    ! mean neutrinosphere energy
REAL(KIND=double), DIMENSION(nnu+1,ny,nz)       :: r_gain          ! gain radius
REAL(KIND=double), DIMENSION(nx,6,ny,nz)        :: r_comp          ! reduced composition mass fraction array
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: dudt_nucMD      ! energy generation rate by nuclear reactions (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: dudt_nuMD       ! energy deposition rate by neutrinos (ergs g^{-1} s^{1})
REAL(KIND=double), DIMENSION(ny,nz)             :: tau_advMD       ! advection time scale
REAL(KIND=double), DIMENSION(ny,nz)             :: tau_heat_nuMD   ! neutrino heating time scale (s)
REAL(KIND=double), DIMENSION(ny,nz)             :: tau_heat_nucMD  ! nuclear heating time scale (s)
REAL(KIND=double), DIMENSION(ny,nz)             :: n_growMD        ! number of convective e-foldings in fluid from shock to gain radius
REAL(KIND=double)                               :: wBV             ! Brunt-Vaisala frequency (positive if unstable) [s^{-1}]
REAL(KIND=double)                               :: twBV            ! Brunt-Vaisala growth time (positive if unstable) [s]
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: wBV_in          ! Brunt-Vaisala frequency (positive if unstable) [s^{-1}]
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: twBV_in         ! Brunt-Vaisala growth time (positive if unstable) [s]
REAL(KIND=double), DIMENSION(nx)                :: wBV_s           ! Brunt-Vaisala frequency due to entropy (positive if unstable) [s^{-1}]
REAL(KIND=double), DIMENSION(nx)                :: twBV_s          ! Brunt-Vaisala growth time due to entropy (positive if unstable) [s]
REAL(KIND=double), DIMENSION(nx)                :: wBV_yl          ! Brunt-Vaisala frequency due to lepton fraction (positive if unstable) [s^{-1}]
REAL(KIND=double), DIMENSION(nx)                :: twBV_yl         ! Brunt-Vaisala growth time due to lepton fraction (positive if unstable) [s]
REAL(KIND=double), DIMENSION(nx)                :: sig             ! set of standard deviations for linear fits
REAL(KIND=double)                               :: a_x             ! parameter a_x in the straight line fit y = a_x x + a_c
REAL(KIND=double)                               :: a_c             ! parameter a_c in the straight line fit y = a_x x + a_c
REAL(KIND=double), DIMENSION(nx,2)              :: rnnu            ! neutrino number density
REAL(KIND=double), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: yl_in           ! lepton fraction
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: wBVMD           ! Brunt Vaisala frequency (positive if unstable) [s^{-1}]
REAL(KIND=double), DIMENSION(nx,ny,nz)          :: ylMD            ! lepton fractionm

REAL(KIND=double), DIMENSION(nx)                :: dr_c            ! radial zone thickness [cm]
REAL(KIND=double), DIMENSION(nx)                :: dvol            ! radial zone volume (cm^{3})
REAL(KIND=double), DIMENSION(nx)                :: m_neut_e        ! unshifted zone_edged Newtonian mass (g)
REAL(KIND=double), DIMENSION(nx)                :: m_neut_c        ! unshifted zone_centered Newtonian mass (g)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: send_buf_nx  ! send buffer for gathering data from processors
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: recv_buf_nx  ! receive buffer for gathering data from processors
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: send_buf_nnu ! send buffer for gathering data from processors
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: recv_buf_nnu ! receive buffer for gathering data from processors

REAL(KIND=double)                               :: x_he            ! mass fraction of helium
REAL(KIND=double)                               :: y_tot           ! abundance fraction of nucleons and nuclei
REAL(KIND=double)                               :: e               ! energy per unit mass (ergs g^{-1})
REAL(KIND=double)                               :: E_env           ! total binding energy of heating region
REAL(KIND=double)                               :: Q_nu            ! total neutrino heating rate in heating region
REAL(KIND=double)                               :: Q_nuc           ! total nuclear heating rate in heating region

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in edit_MD_inout')
 2001 FORMAT (' Deallocation problem for array ',a10,' in edit_MD_inout')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!              \\\\\ LOAD VARIABLES INTO MODULES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer integer scalar parameters
!-----------------------------------------------------------------------

ncycle                    = ncycle_in

!-----------------------------------------------------------------------
!  Transfer real scalar parameters
!-----------------------------------------------------------------------

dtnph                     = dtnph_in
time                      = time_in
t_bounce                  = t_bounce_in
t_tb                      = zero

!-----------------------------------------------------------------------
!
!               \\\\\ DETERMINE WHETHER TO EDIT /////
!               /////   CYCLE NUMBER CRITERIA   \\\\\
!
!-----------------------------------------------------------------------

l_editMDu                 = .false.

IF ( iedMDu /= 0 ) THEN
  nedMDu                  = nedMDu + 1
  IF ( nedMDu >= intedMDu ) THEN
    nedMDu                = 0
    l_editMDu             = .true.
  END IF ! nedMDu >= intedMDu
END IF ! iedMDu /= 0

l_editMDv                 = .false.

IF ( iedMDv /= 0 ) THEN
  nedMDv                  = nedMDv + 1
  IF ( nedMDv >= intedMDv ) THEN
    nedMDv                = 0
    l_editMDv             = .true.
  END IF ! nedMDv >= intedMDv
END IF ! iedMDv /= 0

l_editMDw                 = .false.

IF ( iedMDw /= 0 ) THEN
  nedMDw                  = nedMDw + 1
  IF ( nedMDv >= intedMDv ) THEN
    nedMDw                = 0
    l_editMDw             = .true.
  END IF ! nedMDw >= intedMDw
END IF ! iedMDw /= 0

l_editMDs                 = .false.

IF ( iedMDs /= 0 ) THEN
  nedMDs                  = nedMDs + 1
  IF ( nedMDs >= intedMDs ) THEN
    nedMDs                = 0
    l_editMDs             = .true.
  END IF ! nedMDs >= intedMDs
END IF ! iedMDs /= 0

l_editMDd                 = .false.

IF ( iedMDd /= 0 ) THEN
  nedMDd                  = nedMDd + 1
  IF ( nedMDd >= intedMDd ) THEN
    nedMDd                = 0
    l_editMDd             = .true.
  END IF ! nedMDd >= intedMDd
END IF ! iedMDd /= 0

l_editMDe                 = .false.

IF ( iedMDe /= 0 ) THEN
  nedMDe                  = nedMDe + 1
  IF ( nedMDe >= intedMDe ) THEN
    nedMDe                = 0
    l_editMDe             = .true.
  END IF ! nedMDe >= intedMDe
END IF ! iedMDe /= 0

l_editMDp                 = .false.

IF ( iedMDp /= 0 ) THEN
  nedMDp                  = nedMDp + 1
  IF ( nedMDp >= intedMDp ) THEN
    nedMDp                = 0
    l_editMDp             = .true.
  END IF ! nedMDp >= intedMDp
END IF ! iedMDp /= 0

l_editMDenu               = .false.

IF ( iedMDenu /= 0 ) THEN
  nedMDenu                = nedMDenu + 1
  IF ( nedMDenu >= intedMDenu ) THEN
    nedMDenu              = 0
    l_editMDenu           = .true.
  END IF ! nedMDenu >= intedMDenu
END IF ! iedMDenu /= 0

l_editMDfnu               = .false.

IF ( iedMDfnu /= 0 ) THEN
  nedMDfnu                = nedMDfnu + 1
  IF ( nedMDfnu >= intedMDfnu ) THEN
    nedMDfnu              = 0
    l_editMDfnu           = .true.
  END IF ! nedMDfnu >= intedMDfnu
END IF ! iedMDfnu /= 0

l_editMDa                 = .false.

IF ( iedMDa /= 0 ) THEN
  nedMDa                  = nedMDa + 1
  IF ( nedMDa >= intedMDa ) THEN
    nedMDa                = 0
    l_editMDa             = .true.
  END IF ! nedMDa >= intedMDa
END IF ! iedMDa /= 0

l_editMDx                 = .false.

IF ( iedMDx /= 0 ) THEN
  nedMDx                  = nedMDx + 1
  IF ( nedMDx >= intedMDx ) THEN
    nedMDx                = 0
    l_editMDx             = .true.
  END IF ! nedMDx >= intedMDx
END IF ! iedMDx /= 0

l_editMDye                = .false.

IF ( iedMDye /= 0 ) THEN
  nedMDye                 = nedMDye + 1
  IF ( nedMDye >= intedMDye ) THEN
    nedMDye               = 0
    l_editMDye            = .true.
  END IF ! nedMDye >= intedMDye
END IF ! iedMDye /= 0

l_editMDcm                 = .false.

IF ( iedMDcm /= 0 ) THEN
  nedMDcm                  = nedMDcm + 1
  IF ( nedMDcm >= intedMDcm ) THEN
    nedMDcm                = 0
    l_editMDcm             = .true.
  END IF ! nedMDcm >= intedMDcm
END IF ! iedMDcm /= 0

l_editMDnu                 = .false.

IF ( iedMDnu /= 0 ) THEN
  nedMDnu                  = nedMDnu + 1
  IF ( nedMDnu >= intedMDnu ) THEN
    nedMDnu                = 0
    l_editMDnu             = .true.
  END IF ! nedMDnu >= intedMDnu
END IF ! iedMDnu /= 0

l_editMDnc                 = .false.

IF ( iedMDnc /= 0 ) THEN
  nedMDnc                  = nedMDnc + 1
  IF ( nedMDnc >= intedMDnc ) THEN
    nedMDnc                = 0
    l_editMDnc             = .true.
  END IF ! nedMDnc >= intedMDnc
END IF ! iedMDnc /= 0

l_editMDnl                = .false.

IF ( iedMDnl /= 0 ) THEN
  nedMDnl                 = nedMDnl + 1
  IF ( nedMDnl >= intedMDnl ) THEN
    nedMDnl               = 0
    l_editMDnl            = .true.
  END IF ! nedMDnl >= intedMDnl
END IF ! iedMDnl /= 0

l_editMDne                = .false.

IF ( iedMDne /= 0 ) THEN
  nedMDne                 = nedMDne + 1
  IF ( nedMDne >= intedMDne ) THEN
    nedMDne               = 0
    l_editMDne            = .true.
  END IF ! nedMDne >= intedMDne
END IF ! iedMDne /= 0

l_editMDgx                = .false.

IF ( iedMDgx /= 0 ) THEN
  nedMDgx                 = nedMDgx + 1
  IF ( nedMDgx >= intedMDgx ) THEN
    nedMDgx               = 0
    l_editMDgx            = .true.
  END IF ! nedMDgx >= intedMDgx
END IF ! iedMDgx /= 0

l_editMDgy                = .false.

IF ( iedMDgy /= 0 ) THEN
  nedMDgy                 = nedMDgy + 1
  IF ( nedMDgy >= intedMDgy ) THEN
    nedMDgy               = 0
    l_editMDgy            = .true.
  END IF ! nedMDgy >= intedMDgy
END IF ! iedMDgy /= 0

l_editMDgz                = .false.

IF ( iedMDgz /= 0 ) THEN
  nedMDgz                 = nedMDgz + 1
  IF ( nedMDgz >= intedMDgz ) THEN
    nedMDgz               = 0
    l_editMDgz            = .true.
  END IF ! nedMDgz >= intedMDgz
END IF ! iedMDgz /= 0

!-----------------------------------------------------------------------
!
!               \\\\\ DETERMINE WHETHER TO EDIT /////
!               /////       TIME CRITERIA       \\\\\
!
!-----------------------------------------------------------------------

l_MDedit                  = .false.

IF ( i_editMD /= 0 ) THEN

  IF ( t_bounce <= zero ) THEN

!-----------------------------------------------------------------------
!  Initialize time criterion
!-----------------------------------------------------------------------

    IF ( first ) THEN
      tmult               = 1.d+3/dt_MDedit1
      itimeprev           = int( time * tmult )
      first               = .false.
    END IF

!-----------------------------------------------------------------------
!  Generate a MD edit every multiple of dt_MDedit1 ms up to bounce bounce
!-----------------------------------------------------------------------

    tmult                 = 1.d+3/dt_MDedit1
    itime                 = int( time * tmult )
    IF ( itime /= itimeprev ) THEN
      l_MDedit            = .true.
      itimeprev           = itime
    END IF
    
  ELSE

    t_tb                  = time - t_bounce

!-----------------------------------------------------------------------
!  Generate an MD edit at bounce
!-----------------------------------------------------------------------

    IF ( t_tb - dtnph <= zero  .and.  t_tb > zero ) l_MDedit = .true.

!-----------------------------------------------------------------------
!  Generate an MD edit 1 ms after bounce
!-----------------------------------------------------------------------

    IF ( t_tb - dtnph <= 1.d-3  .and.  t_tb > 1.d-3 ) l_MDedit = .true.

!-----------------------------------------------------------------------
!  Generate an MD edit 10 ms after bounce
!-----------------------------------------------------------------------

    IF ( t_tb - dtnph <= 1.d-2  .and.  t_tb > 1.d-2 ) l_MDedit = .true.

!-----------------------------------------------------------------------
!  Initialize time criterion
!-----------------------------------------------------------------------

    IF ( first ) THEN
      tmult               = 1.d+3/dt_MDedit2
      itimeprev           = int( t_tb * tmult )
      first               = .false.
    END IF

!-----------------------------------------------------------------------
!  Generate an MD edit every multiple of dt_MDedit2 ms after bounce
!-----------------------------------------------------------------------

    tmult                 = 1.d+3/dt_MDedit2
    itime                 = int( t_tb * tmult )
    IF ( itime /= itimeprev ) THEN
      l_MDedit            = .true.
      itimeprev           = itime
    END IF

  END IF

END IF ! i_editMD /= 0

!-----------------------------------------------------------------------
!  Return if l_MDedit = .false
!-----------------------------------------------------------------------

IF (       ( .not. l_editMDu  ) &
&   .and.  ( .not. l_editMDv  ) &
&   .and.  ( .not. l_editMDw  ) &
&   .and.  ( .not. l_editMDs  ) &
&   .and.  ( .not. l_editMDd  ) &
&   .and.  ( .not. l_editMDe  ) &
&   .and.  ( .not. l_editMDp  ) &
&   .and.  ( .not. l_editMDenu) &
&   .and.  ( .not. l_editMDfnu) &
&   .and.  ( .not. l_editMDa  ) &
&   .and.  ( .not. l_editMDx  ) &
&   .and.  ( .not. l_editMDye ) &
&   .and.  ( .not. l_editMDcm ) &
&   .and.  ( .not. l_editMDnu ) &
&   .and.  ( .not. l_editMDnc ) &
&   .and.  ( .not. l_editMDnl ) &
&   .and.  ( .not. l_editMDne ) &
&   .and.  ( .not. l_editMDgx ) &
&   .and.  ( .not. l_editMDgy ) &
&   .and.  ( .not. l_editMDgz ) &
&   .and.  ( .not. l_MDedit   ) ) RETURN

!-----------------------------------------------------------------------
!
!                \\\\\ COMPUTE VARIABLES TO EDIT /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  MGFLD shifted indices
!-----------------------------------------------------------------------

jr_min                    = imin + 1
jr_max                    = imax + 1
jr_maxp                   = jr_max + 1

!-----------------------------------------------------------------------
!  Mean enclosed mass
!-----------------------------------------------------------------------

dr_c(1:imax)              = x_e_in(2:imax+1) - x_e_in(1:imax)
dvol(1:imax)              = frpi * dr_c(1:imax) * ( x_e_in(1:imax) * ( x_e_in(1:imax) + dr_c(1:imax) ) &
&                         + dr_c(1:imax) * dr_c(1:imax) * third ) 

i                         = imin
IF ( x_e_in(i) == zero ) THEN
  m_neut_e(i)             = zero
ELSE
  m_neut_e(i)             = rhobar(i) * frpith * x_e_in(i)**3
END IF

DO  i = imin,imax
  m_neut_e(i+1)           = m_neut_e(i) + rhobar(i) * dvol(i)
END  DO

m_neut_c(imin:imax)       = half * ( m_neut_e(imin:imax ) + m_neut_e(imin+1:imax+1) )

!-----------------------------------------------------------------------
!  Shock location
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    CALL findshock( imin+1, imax, j, k, pqmin, j_shock, j_shock_mx, j_shock_mn, &
&    nx, x_c_in, m_neut_c, r_shk, m_shk )
    i_shock                 = MAX( j_shk_radial_p(j,k) - 1, 1 )
    i_shock_in(j,k)         = i_shock
    i_shock_mx              = MAX( j_shock_mx - 1, 1 )
    i_shock_mn              = MAX( j_shock_mn - 1, 1 )
    j_shock_in(j,k)         = i_shock + 1
    r_shock_in1             = x_c_in(i_shock)
    r_shock_in_mn(j,k)      = x_c_in(i_shock_mn)
    r_shock_in_mx(j,k)      = x_c_in(i_shock_mx)
    IF ( i_shock_mn == 1 ) r_shock_in_mn(j,k) = 0.d0
    IF ( i_shock_mx == 1 ) r_shock_in_mx(j,k) = 0.d0
    IF ( i_shock    == 1 ) THEN
      r_shock_in (j,k)      = 0.d0
    ELSE
      IF ( ( pq_x(j_shock+1,j,k) + pqy_x(j_shock+1,j,k) )/( aesv(j_shock+1,1,j,k) + epsilon )  &
&        > ( pq_x(j_shock-1,j,k) + pqy_x(j_shock-1,j,k) )/( aesv(j_shock-1,1,j,k) + epsilon ) ) THEN
        i_shockm            = i_shock + 1
      ELSE
        i_shockm            = i_shock - 1
      END IF ! pq_x/aesv
      r_shock_in2           = x_c_in(i_shockm)
      q_shock_in1           = ( pq_x(i_shock+1 ,j,k) + pqy_x(i_shock+1 ,j,k) )/( aesv(i_shock+1 ,1,j,k) + epsilon )
      q_shock_in2           = ( pq_x(i_shockm+1,j,k) + pqy_x(i_shockm+1,j,k) )/( aesv(i_shockm+1,1,j,k) + epsilon )
      r_shock_in(j,k)       = mean( q_shock_in1, q_shock_in2, r_shock_in1, r_shock_in2 )
    END IF ! i_shock   -1 /= 1
  END DO ! j = 1,ij_ray_dim
END DO ! k = 1,ik_ray_dim

!-----------------------------------------------------------------------
!  Neutrino luminosities
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    CALL luminosity_MD( imin, imax, j, k, ij_ray_dim, ik_ray_dim, nx, nez, &
&    nnu, x_e_in, psi1_in, unue_in, dunue_in, lum_in )
  END DO ! j = 1,ij_ray_dim
END DO ! k = 1,ik_ray_dim

!-----------------------------------------------------------------------
!  Neutrino rms energies
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    CALL e_rms_MD( imin, imax, j, k, ij_ray_dim, ik_ray_dim, nx, nez, nnu, &
&    psi0_in, psi1_in, unu_in, dunu_in, e_rms_stat_in, e_rms_trns_in )
  END DO ! j = 1,ij_ray_dim
END DO ! k = 1,ik_ray_dim

!-----------------------------------------------------------------------
!  Mean neutrinospheres
!-----------------------------------------------------------------------

rstmss(1)                 = zero
DO i = imin+1,imax+1
  rstmss(i)               = rstmss(i-1) + frpith * ( x_e_in(i)**3 - x_e_in(i-1)**3 ) * rhobar(i-1)
END DO

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim

    rho(jr_min:jr_max)    = rhoMD_in(imin:imax,j,k)
    t  (jr_min:jr_max)    = tMD_in  (imin:imax,j,k)

    CALL nu_sphere( jr_min, jr_maxp, j, k, ij_ray_dim, ik_ray_dim, x_e_in, &
&    rho, t, rstmss, nx, nez, nnu, j_sphere, r_sphere, d_sphere, t_sphere, &
&    m_sphere )

    CALL nu_sphere_mean( jr_min, jr_maxp, j, k, ij_ray_dim, ik_ray_dim,    &
&    e_rms_trns_in, j_sphere, r_sphere, d_sphere, t_sphere, m_sphere, nx,  &
&    nez, nnu, rsphere_mean_in, dsphere_mean_in, tsphere_mean_in,          &
&    msphere_mean_in, esphere_mean_in, jsphere_mean_in )

  END DO ! j = 1,ij_ray_dim
END DO ! k = 1,ik_ray_dim

!-----------------------------------------------------------------------
!  Gain radius
!
!   Determine r_gain(1:2,i_ray) by working inward from the shock and
!    locating the radius where dunujeadt changes sign
!
!   Determine r_gain(5,i_ray)  by working inward from the shock and
!    locating the radius where two adjecent dudt_nu are positive
!    followed inward by two adjecent negative dudt_nu
!
!   Set gain radius to zero if either
!    (1) neutrinosphere == 0
!    (2) no shock
!-----------------------------------------------------------------------

DO n = 1,2
  DO j = 1,ij_ray_dim
    DO k = 1,ik_ray_dim
      j_gain              = jr_min
      i_gain              = jr_min - 1
      IF ( jsphere_mean_in(n,j,k) == jr_min  .or.  r_shock_in(j,k) == zero ) THEN
        r_gain_in(n,j,k)  = zero
      ELSE
        DO i = jr_max,jsphere_mean_in(n,j,k)+1,-1
          IF ( nse_in(i+1,j,k) == 0 ) CYCLE
          IF ( dunujeadt(i,n,j,k) > zero  .and.  dunujeadt(i+1,n,j,k) < zero ) THEN
            j_gain        = i
            i_gain        = i - 1
            EXIT
          END IF ! dunujeadt(i,n,j,k) > zero  .and.  dunujeadt(i+1,n,j,k) < zero
        END DO ! i = jr_max,jsphere_mean_in(n,j,k)+1,-1
        IF ( j_gain == jr_min ) THEN
          r_gain_in(n,j,k)  = zero
          CYCLE
        END IF ! j_gain == jr_min
        r_gain_in(n,j,k)  = mean( dunujeadt(j_gain,n,j,k), -dunujeadt(j_gain+1,n,j,k), &
&                                 x_c_in(i_gain+1), x_c_in(i_gain) )
      END IF ! jsphere_mean_in(n,j) == jr_min  .or.  r_shock_in(j) == zero
    END DO ! k = 1,ik_ray_dim
  END DO ! j = 1,ij_ray_dim
END DO ! n = 1,2

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    j_gain                = jr_min
    i_gain                = jr_min - 1
    i_gain_in(j,k)        = i_gain
    IF ( jsphere_mean_in(1,j,k) == jr_min  .or.                        &
&        jsphere_mean_in(2,j,k) == jr_min  .or.                        &
&        r_shock_in(j,k) == zero )                                      THEN
      r_gain_in(5,j,k)    = zero
    ELSE
      DO i = jr_max,jsphere_mean_in(2,j,k)+2,-1
        IF ( nse_in(I+1,j,k) == 0 ) CYCLE
        IF ( dudt_nu(i-1,j,k) < zero  .and.                            &
&            dudt_nu(i,j,k) < zero    .and.                            &
&            dudt_nu(i+1,j,k) > zero  .and.                            &
&            dudt_nu(i+2,j,k) > zero )                                  THEN
          j_gain          = i
          i_gain          = i - 1
          i_gain_in(j,k)  = i_gain
          EXIT
        END IF ! dudt_nu(i,j,k) > zero  .and.  dudt_nu(i+1,j,k) < zero 
      END DO ! i = jr_max,jsphere_mean_in(2,j,k)+1,-1
      IF ( j_gain == jr_min ) THEN
        r_gain_in(5,j,k)  = zero
        CYCLE
      ELSE
        r_gain_in(5,j,k)  = mean( -dudt_nu(j_gain,j,k), dudt_nu(j_gain+1,j,k), &
&                                  x_c_in(i_gain+1), x_c_in(i_gain) )
      END IF ! j_gain == jr_min
    END IF ! jsphere_mean_in(1,j,k) == jr_min  .or.  jsphere_mean_in(2,j,k) == jr_min  .or.  r_shock_in(j,k) == zero
  END DO ! k = 1,ik_ray_dim
END DO ! j = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Advection time scale
!
!   Determine the advection time scale by integrating dr/(-u(r)) from
!    the neutrinosphere to the shock
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    IF ( jsphere_mean_in(1,j,k) == jr_min              .or.             &
&        jsphere_mean_in(2,j,k) == jr_min              .or.             &
&        r_shock_in(j,k)        == zero                .or.             &
&        i_shock_in(j,k) - 2    <= i_gain_in(j,k) + 2  .or.             &
&        i_gain_in(j,k)         == 1 )                   THEN
      tau_adv_in(j,k)     = 1.d+100
    ELSE
      nwt                 = 0
      sig                 = 1.d0
      CALL linear_fit( x_c_in, uMD_in(:,j,k), nx, i_gain_in(j,k), i_shock_in(j,k)-2, sig, nwt, a_c, a_x )
      tau_adv_in(j,k)     = zero
      DO i = i_shock_in(j,k) - 2, i_gain_in(j,k) + 1, -1
        tau_adv_in(j,k)   = tau_adv_in(j,k) + dr_c(i)/( DABS( a_x * x_c_in(i) + a_c ) + epsilon )
      END DO ! i = i_shock_in(j,k),i_gain_in(j,k),-1
    END IF ! jsphere_mean_in(1,j,k) == jr_min. etc.
  END DO ! k = 1,ik_ray_dim
END DO ! j = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Neutrino heating time scale
!
!   Determine the heating time scale by integrating E/dEdt from
!    the neutrinosphere to the shock
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    IF ( tau_adv_in(j,k) == 1.d+100 ) THEN
      tau_heat_nu_in(j,k) = 1.d+100
      tau_heat_nuc_in(j,k) = 1.d+100
    ELSE
      E_env               = zero
      Q_nu                = zero
      Q_nuc               = zero
      DO i = i_shock_in(j,k) - 2, i_gain_in(j,k) + 1, -1
        e                 = half * ( uMD_in(i,j,k) * uMD_in(i,j,k)      &
&                         + vMD_in(i,j,k) * vMD_in(i,j,k) )             &
&                         + aesv(i+1,2,j,k) - e_bar0 + gtot_pot_c(i,j,k)
        E_env             = E_env + e * ( m_neut_e(i+1) - m_neut_e(i) )
        Q_nu              = Q_nu  + dudt_nu(i+1,j,k)  * ( m_neut_e(i+1) - m_neut_e(i) )
        Q_nuc             = Q_nuc + dudt_nuc(i+1,j,k) * ( m_neut_e(i+1) - m_neut_e(i) )
      END DO ! i = i_shock_in(j,k),i_gain_in(j,k)+1,-1
      tau_heat_nu_in(j,k) = E_env/( Q_nu + epsilon )
      tau_heat_nuc_in(j,k) = E_env/( Q_nuc + epsilon )
    END IF ! tau_adv_in(j,k) == 1.d+100
  END DO ! k = 1,ik_ray_dim
END DO ! j = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Convective growth time scale
!
!   Determine the convective growth time scale by integrating E/dEdt from
!    the neutrinosphere to the shock
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    IF ( tau_adv_in(j,k) == 1.d+100 ) THEN
      n_grow_in(j,k)      = 1.d+100
    ELSE
      nwt                 = 0
      sig                 = 1.d0
      CALL linear_fit( x_c_in, uMD_in(:,j,k), nx, i_gain_in(j,k), i_shock_in(j,k)-2, sig, nwt, a_c, a_x )
      n_grow_in(j,k)      = zero
      DO i = i_shock_in(j,k) - 2, i_gain_in(j,k) + 1, -1
        CALL Brunt_Vaisala_ye( i, j, k, nx, ij_ray_dim, ik_ray_dim, x_c_in, rhoMD_in, yeMD_in, &
&        grav_x_c_in, wBV, twBV )
        n_grow_in(j,k)    = n_grow_in(j,k) + DMAX1( wBV, zero ) * dr_c(i)/( DABS( a_x * x_c_in(i) + a_c ) + epsilon )
      END DO ! i = i_shock_in(j) - 2, i_gain_in(j) + 1, -1
    END IF ! tau_adv_in(j) == 1.d+100
  END DO ! k = 1,ik_ray_dim
END DO ! j = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Sound velocity
!-----------------------------------------------------------------------

DO k = 1,ik_ray_dim
  DO j = 1,ij_ray_dim
    DO i = imin,imax
      v_csoundMD_2        = ( uMD_in(i,j,k) * uMD_in(i,j,k)          &
&                         +   vMD_in(i,j,k) * vMD_in(i,j,k)          &
&                         +   wMD_in(i,j,k) * wMD_in(i,j,k) + epsilon ) &
&                         * rhoMD_in(i,j,k)/( aesv(i+1,12,j,k) * aesv(i+1,1,j,k) )
      v_csoundMD_in(i,j,k) = DSQRT( v_csoundMD_2 )
    END DO ! i = imin,imax
  END DO ! j = 1,ij_ray_dim
END DO ! k = 1,ik_ray_dim

!-----------------------------------------------------------------------
!  Neutrino energz density
!-----------------------------------------------------------------------

e_nud_MD_in(imin:imax,:,:) = e_nu_MD_in(imin:imax,:,:)/rhoMD_in(imin:imax,:,:)

!-----------------------------------------------------------------------
!  Mean nuclear mass number
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    DO i = imin,imax
      IF ( nse_in(i,j,k) == 1 ) THEN
        x_he              = DMAX1( 1.d0 - aesv(i+1,7,j,k) - aesv(i+1,8,j,k) - aesv(i+1,9,j,k), zero )
        y_tot             = aesv(i+1,7,j,k) + aesv(i+1,8,j,k) + x_he/4.d0 + aesv(i+1,9,j,k)/( aesv(i+1,10,j,k) + epsilon )
        aMD_in(i,j,k)     = 1.d0/( y_tot + epsilon )
      ELSE
        y_tot             = SUM( xn_in(i,1:nnc-1,j,k)/( a_nuc(1:nnc-1) + epsilon ) )
        y_tot             = y_tot + xn_in(i,nnc,j,k)/( a_nuc_rep(i+1) + epsilon ) 
        aMD_in(i,j,k)     = 1.d0/( y_tot + epsilon )
      END IF
    END DO ! i = imin,imax
  END DO ! k = kmin,kmax
END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!  Ledoux stability
!
!   Determine the Ledoux stability, coupling neutrinos to matter at high
!    density, uncoupling them at low density
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim

    rho(jr_min:jr_max)    = rhoMD_in(imin:imax,j,k)
    t  (jr_min:jr_max)    = tMD_in  (imin:imax,j,k)
    ye (jr_min:jr_max)    = yeMD_in (imin:imax,j,k)
    r  (imin:imax+1)      = x_c_in  (imin:imax+1)
    g_acc_e(imin:imax+1)  = grav_x_e_in(imin:imax+1,j,k)
  

    wBV_in                = zero
    twBV_in               = zero
    wBV_s                 = zero
    twBV_s                = zero
    wBV_yl                = zero
    twBV_yl               = zero

    CALL esrgneldnu( jr_min, jr_max, nx, rho, t, ye, j, k, ye_yl )
    CALL eqstld( jr_min, jr_max, rho, t, ye_yl, j, k )

    DO jr = jr_min, jr_max
      CALL Brunt_Vaisala_yl( jr, j, k, nx, r, rho, ye_yl, g_acc_e, wBV_in(jr,j,k), &
&      twBV_in(jr,j,k), wBV_s(jr), twBV_s(jr), wBV_yl(jr), twBV_yl(jr) )
    END DO ! jr = jr_min, jr_max

  END DO ! k = 1,ik_ray_dim
END DO ! j = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Lepton fraction 
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    rnnu                  = zero
    DO i = imin, imax
      DO n = 1, 2
        IF ( nnugp(n) /= 0 ) THEN
          rnnu(i,n)       = SUM( unu_in(i,:,j,k) * unu_in(i,:,j,k) * dunu_in(i,:,j,k) * psi0_in(i,:,n,j,k) ) * ( ncoef/rhoMD_in(i,j,k) )
        END IF !  nnugp(n) ne 0
      END DO !  n = 1, 2
      yl_in(i,j,k)        = yeMD_in(i,j,k) + ( rnnu(i,1) - rnnu(i,2) ) * rmu
    END DO ! i = imin, imax
  END DO ! k = 1,ik_ray_dim
END DO ! j = 1,ij_ray_dim

!-----------------------------------------------------------------------
!
!       \\\\\ GATHER VARIABLES TO EDIT FROM ALL PROCESSORS /////
!
!-----------------------------------------------------------------------

uMD                       = zero
vMD                       = zero
wMD                       = zero
rhoMD                     = zero
sMD                       = zero
yeMD                      = zero
v_csoundMD                = zero
e_nu_MD                   = zero
f_nu_MD                   = zero
grav_x_cMD                = zero
grav_y_cMD                = zero
grav_z_cMD                = zero
lumMD                     = zero
e_rmsMD                   = zero
eMD                       = zero
pMD                       = zero
aMD                       = zero
dudt_nuc                  = zero
dudt_nu                   = zero
aesvMD                    = zero
wBVMD                     = zero
ylMD                      = zero

!-----------------------------------------------------------------------
!  Allocate send and receive buffer
!-----------------------------------------------------------------------

ALLOCATE (send_buf_nx(16,nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'send_buf_nnu  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (recv_buf_nx(16,nx,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'recv_buf_nnu  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Load send buffer
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    send_buf_nx( 1,1:nx,j,k) = uMD_in       (1:nx,j,k)
    send_buf_nx( 2,1:nx,j,k) = vMD_in       (1:nx,j,k)
    send_buf_nx( 3,1:nx,j,k) = wMD_in       (1:nx,j,k)
    send_buf_nx( 4,1:nx,j,k) = rhoMD_in     (1:nx,j,k)
    send_buf_nx( 5,1:nx,j,k) = yeMD_in      (1:nx,j,k)
    send_buf_nx( 6,1:nx,j,k) = v_csoundMD_in(1:nx,j,k)
    send_buf_nx( 7,1:nx,j,k) = e_nud_MD_in  (1:nx,j,k)
    send_buf_nx( 8,1:nx,j,k) = f_nu_MD_in   (1:nx,j,k)
    send_buf_nx( 9,1:nx,j,k) = grav_x_c_in  (1:nx,j,k)
    send_buf_nx(10,1:nx,j,k) = grav_y_c_in  (1:nx,j,k)
    send_buf_nx(11,1:nx,j,k) = grav_z_c_in  (1:nx,j,k)
    send_buf_nx(12,1:nx,j,k) = aMD_in       (1:nx,j,k)
    send_buf_nx(13,1:nx,j,k) = dudt_nuc     (1:nx,j,k)
    send_buf_nx(14,1:nx,j,k) = dudt_nu      (1:nx,j,k)
    send_buf_nx(15,1:nx,j,k) = wBV_in       (1:nx,j,k)
    send_buf_nx(16,1:nx,j,k) = yl_in        (1:nx,j,k)
  END DO ! k = 1,ik_ray_dim
END DO ! j = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Gather quantities
!-----------------------------------------------------------------------

c_gath_recv               = 16 * nx * ij_ray_dim * ik_ray_dim
c_gath_send               = 16 * nx * ij_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( send_buf_nx, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf_nx, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!                ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Unpack receive buffer
!-----------------------------------------------------------------------

  DO m = 0,n_proc-1
    mj                    = MOD( m, n_proc_y )
    mk                    = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                 = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk               = mj * ij_ray_dim + j
        uMD       (1:nx,jsk,ksk) = recv_buf_nx( 1,1:nx,j,k,m+1)
        vMD       (1:nx,jsk,ksk) = recv_buf_nx( 2,1:nx,j,k,m+1)
        wMD       (1:nx,jsk,ksk) = recv_buf_nx( 3,1:nx,j,k,m+1)
        rhoMD     (1:nx,jsk,ksk) = recv_buf_nx( 4,1:nx,j,k,m+1)
        yeMD      (1:nx,jsk,ksk) = recv_buf_nx( 5,1:nx,j,k,m+1)
        v_csoundMD(1:nx,jsk,ksk) = recv_buf_nx( 6,1:nx,j,k,m+1)
        e_nu_MD   (1:nx,jsk,ksk) = recv_buf_nx( 7,1:nx,j,k,m+1)
        f_nu_MD   (1:nx,jsk,ksk) = recv_buf_nx( 8,1:nx,j,k,m+1)
        grav_x_cMD(1:nx,jsk,ksk) = recv_buf_nx( 9,1:nx,j,k,m+1)
        grav_y_cMD(1:nx,jsk,ksk) = recv_buf_nx(10,1:nx,j,k,m+1)
        grav_z_cMD(1:nx,jsk,ksk) = recv_buf_nx(11,1:nx,j,k,m+1)
        aMD       (1:nx,jsk,ksk) = recv_buf_nx(12,1:nx,j,k,m+1)
        dudt_nucMD(1:nx,jsk,ksk) = recv_buf_nx(13,1:nx,j,k,m+1)
        dudt_nuMD (1:nx,jsk,ksk) = recv_buf_nx(14,1:nx,j,k,m+1)
        wBVMD     (1:nx,jsk,ksk) = recv_buf_nx(15,1:nx,j,k,m+1)
        ylMD      (1:nx,jsk,ksk) = recv_buf_nx(16,1:nx,j,k,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!                 ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Deallocate send and receive buffer
!-----------------------------------------------------------------------

DEALLOCATE (send_buf_nx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'send_buf_nx   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (recv_buf_nx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'recv_buf_nx   '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate send and receive buffer
!-----------------------------------------------------------------------

ALLOCATE (send_buf_nx(nx,nnc,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'send_buf_nnu  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (recv_buf_nx(nx,nnc,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'recv_buf_nnu  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Load send buffer
!-----------------------------------------------------------------------

DO ic = 1, nnc
  DO j = 1,ij_ray_dim
    DO k = 1,ik_ray_dim
      send_buf_nx(1:nx,ic,j,k) = xn_in(1:nx,ic,j,k)
    END DO ! k = 1,ik_ray_dim
  END DO ! j = 1,ij_ray_dim
END DO ! ic = 1, 12

!-----------------------------------------------------------------------
!  Gather quantities
!-----------------------------------------------------------------------

c_gath_recv               = nx * nnc * ij_ray_dim * ik_ray_dim
c_gath_send               = nx * nnc * ij_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( send_buf_nx, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf_nx, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!                ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Unpack receive buffer
!-----------------------------------------------------------------------

  DO m = 0,n_proc-1
    mj                    = MOD( m, n_proc_y )
    mk                    = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                 = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk               = mj * ij_ray_dim + j
        DO ic = 1, nnc
          xnMD(1:nx,ic,jsk,ksk) = recv_buf_nx(1:nx,ic,j,k,m+1)
        END DO ! ic = 1, nnc
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!                 ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Deallocate send and receive buffer
!-----------------------------------------------------------------------

DEALLOCATE (send_buf_nx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'send_buf_nx   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (recv_buf_nx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'recv_buf_nx   '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate send and receive buffer
!-----------------------------------------------------------------------

ALLOCATE (send_buf_nx(nx,2*nnu+12,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'send_buf_nx  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (recv_buf_nx(nx,2*nnu+12,ij_ray_dim,ik_ray_dim,n_proc_y*n_proc_z), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'recv_buf_nx  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Load send buffer
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    DO ic = 1, nnu
      send_buf_nx(1:nx,ic      ,j,k) = lum_in(1:nx,ic,j,k)
    END DO ! ic = 1, nnu
    DO ic = 1, nnu
      send_buf_nx(1:nx,ic  +nnu,j,k) = e_rms_stat_in(1:nx,ic,j,k)
    END DO ! ic = 1, nnu
    DO ic = 1, 12
      send_buf_nx(1:nx,ic+2*nnu,j,k) = aesv(1:nx,ic,j,k)
    END DO ! ic = 1, 12
  END DO ! k = 1,ik_ray_dim
END DO ! j = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Gather quantities
!-----------------------------------------------------------------------

c_gath_recv               = nx * ( 2 * nnu + 12 ) * ij_ray_dim * ik_ray_dim
c_gath_send               = nx * ( 2 * nnu + 12 ) * ij_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( send_buf_nx, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf_nx, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!                ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Unpack receive buffer
!-----------------------------------------------------------------------

  DO m = 0,n_proc-1
    mj                    = MOD( m, n_proc_y )
    mk                    = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                 = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk               = mj * ij_ray_dim + j
        DO ic = 1, nnu
          lumMD(1:nx,ic,jsk,ksk)   = recv_buf_nx(1:nx,ic      ,j,k,m+1)
        END DO ! ic = 1, nnu
        DO ic = 1, nnu
          e_rmsMD(1:nx,ic,jsk,ksk) = recv_buf_nx(1:nx,ic  +nnu,j,k,m+1)
        END DO ! ic = 1, nnu
        DO ic = 1, 12
          aesvMD(1:nx,ic,jsk,ksk)  = recv_buf_nx(1:nx,ic+2*nnu,j,k,m+1)
        END DO ! ic = 1, 12
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!  Internal energz
!-----------------------------------------------------------------------

  eMD(imin:imax,:,:)             = aesvMD(imin+1:imax+1,2,:,:)

!-----------------------------------------------------------------------
!  Pressure
!-----------------------------------------------------------------------

  pMD(imin:imax,:,:)             = aesvMD(imin+1:imax+1,1,:,:)

!-----------------------------------------------------------------------
!  Entropy
!-----------------------------------------------------------------------

  sMD(imin:imax,:,:)             = aesvMD(imin+1:imax+1,3,:,:)

!-----------------------------------------------------------------------
!                 ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Deallocate send and receive buffer
!-----------------------------------------------------------------------

DEALLOCATE (send_buf_nx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'send_buf_nx   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (recv_buf_nx, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'recv_buf_nx   '; WRITE (nlog,2001) var_name; END IF

!........NSE flag

!-----------------------------------------------------------------------
!  Allocate receive buffer
!-----------------------------------------------------------------------

ALLOCATE (irecv_buf(nx,ij_ray_dim,ij_ray_dim,n_proc_y*n_proc_z), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'irecv_buf     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Gather quantities
!-----------------------------------------------------------------------

c_gath_recv               = nx * ij_ray_dim * ik_ray_dim
c_gath_send               = nx * ij_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( nse_in, c_gath_send, MPI_INTEGER, irecv_buf, c_gath_recv, &
& MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!      ||||| The following is performed on processor 0 only |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

  DO m = 0,n_proc-1
    mj                    = MOD( m, n_proc_y )
    mk                    = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                 = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk               = mj * ij_ray_dim + j
        nseMD(1:nx,jsk,ksk) = irecv_buf(1:nx,j,k,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!                 ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Deallocate send and receive buffer
!-----------------------------------------------------------------------

DEALLOCATE (irecv_buf, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'irecv_buf     '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Allocate send and receive buffer
!-----------------------------------------------------------------------

ALLOCATE (send_buf_nnu((7+6*nnu+1),ij_ray_dim,ij_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'send_buf_nnu  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (recv_buf_nnu((7+6*nnu+1),ij_ray_dim,ij_ray_dim,n_proc_y*n_proc_z), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'recv_buf_nnu  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Load send buffer
!-----------------------------------------------------------------------

DO j = 1,ij_ray_dim
  DO k = 1,ik_ray_dim
    send_buf_nnu(1                  ,j,k) = r_shock_in(j,k)
    send_buf_nnu(2                  ,j,k) = r_shock_in_mn(j,k)
    send_buf_nnu(3                  ,j,k) = r_shock_in_mx(j,k)
    send_buf_nnu(4                  ,j,k) = tau_adv_in(j,k)
    send_buf_nnu(5                  ,j,k) = tau_heat_nu_in(j,k)
    send_buf_nnu(6                  ,j,k) = tau_heat_nuc_in(j,k)
    send_buf_nnu(7                  ,j,k) = n_grow_in(j,k)
    send_buf_nnu(7      +1:7+1*nnu  ,j,k) = rsphere_mean_in(1:nnu,j,k)
    send_buf_nnu(7+1*nnu+1:7+2*nnu  ,j,k) = dsphere_mean_in(1:nnu,j,k)
    send_buf_nnu(7+2*nnu+1:7+3*nnu  ,j,k) = tsphere_mean_in(1:nnu,j,k)
    send_buf_nnu(7+3*nnu+1:7+4*nnu  ,j,k) = msphere_mean(1:nnu,j,k)
    send_buf_nnu(7+4*nnu+1:7+5*nnu  ,j,k) = esphere_mean(1:nnu,j,k)
    send_buf_nnu(7+5*nnu+1:7+6*nnu+1,j,k) = r_gain_in(1:nnu+1,j,k)
  END DO ! ik_ray_dim
END DO ! ij_ray_dim

!-----------------------------------------------------------------------
!  Gather quantities
!-----------------------------------------------------------------------

c_gath_recv               = ( 7 + 6 * nnu + 1 ) * ij_ray_dim * ik_ray_dim
c_gath_send               = ( 7 + 6 * nnu + 1 ) * ij_ray_dim * ik_ray_dim

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER( send_buf_nnu, c_gath_send, MPI_DOUBLE_PRECISION, recv_buf_nnu, c_gath_recv, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------------
!                ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!  Unpack receive buffer
!-----------------------------------------------------------------------

  DO m = 0,n_proc-1
    mj                    = MOD( m, n_proc_y )
    mk                    = m / n_proc_y     
    DO k = 1,ik_ray_dim
      ksk                 = mk * ik_ray_dim + k
      DO j = 1,ij_ray_dim
        jsk               = mj * ij_ray_dim + j
        r_shock       (      jsk,ksk) = recv_buf_nnu(1                  ,j,k,m+1)
        r_shock_mn    (      jsk,ksk) = recv_buf_nnu(2                  ,j,k,m+1)
        r_shock_mx    (      jsk,ksk) = recv_buf_nnu(3                  ,j,k,m+1)
        tau_advMD     (      jsk,ksk) = recv_buf_nnu(4                  ,j,k,m+1)
        tau_heat_nuMD (      jsk,ksk) = recv_buf_nnu(5                  ,j,k,m+1)
        tau_heat_nucMD(      jsk,ksk) = recv_buf_nnu(6                  ,j,k,m+1)
        n_growMD      (      jsk,ksk) = recv_buf_nnu(7                  ,j,k,m+1)
        rsphere_mean(1:nnu  ,jsk,ksk) = recv_buf_nnu(7      +1:7+1*nnu  ,j,k,m+1)
        dsphere_mean(1:nnu  ,jsk,ksk) = recv_buf_nnu(7+1*nnu+1:7+2*nnu  ,j,k,m+1)
        tsphere_mean(1:nnu  ,jsk,ksk) = recv_buf_nnu(7+2*nnu+1:7+3*nnu  ,j,k,m+1)
        msphere_mean(1:nnu  ,jsk,ksk) = recv_buf_nnu(7+3*nnu+1:7+4*nnu  ,j,k,m+1)
        esphere_mean(1:nnu  ,jsk,ksk) = recv_buf_nnu(7+4*nnu+1:7+5*nnu  ,j,k,m+1)
        r_gain      (1:nnu+1,jsk,ksk) = recv_buf_nnu(7+5*nnu+1:7+6*nnu+1,j,k,m+1)
      END DO ! j = 1,ij_ray_dim
    END DO ! k = 1,ik_ray_dim
  END DO ! m = 0,n_proc-1

!-----------------------------------------------------------------------
!                 ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

!-----------------------------------------------------------------------
!  Deallocate send and receive buffer
!-----------------------------------------------------------------------

DEALLOCATE (send_buf_nnu, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'send_buf_nnu  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (recv_buf_nnu, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'recv_buf_nnu  '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!                ||||| MYID == 0 ONLY STARTS HERE |||||
!-----------------------------------------------------------------------

IF ( myid == 0 ) THEN

!-----------------------------------------------------------------------
!
!                      \\\\\ COMPOSITION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Composition indexing
!-----------------------------------------------------------------------

  IF ( first_c ) THEN
    first_c               = .false.
    i_n                   = nuc_number + 1
    i_p                   = nuc_number + 1
    i_4He                 = nuc_number + 1
    i_12C                 = nuc_number + 1
    i_16O                 = nuc_number + 1
    i_20Ne                = nuc_number + 1
    i_24Mg                = nuc_number + 1
    i_28Si                = nuc_number + 1
    i_32S                 = nuc_number + 1
    i_36Ar                = nuc_number + 1
    i_40Ca                = nuc_number + 1
    i_44Ti                = nuc_number + 1
    i_48Cr                = nuc_number + 1
    i_52Fe                = nuc_number + 1
    i_56Ni                = nuc_number + 1
    i_60Zn                = nuc_number + 1
    DO i = 1,nuc_number
      IF (      a_name(i) == ' 12C ' ) THEN
        i_12C             = i
      END IF 
      IF (      a_name(i) == ' 16O ' ) THEN
        i_16O             = i
      END IF 
      IF (      a_name(i) == ' 20Ne' ) THEN
        i_20Ne            = i
      END IF 
      IF (      a_name(i) == ' 24Mg' ) THEN
        i_24Mg            = i
      END IF 
      IF (      a_name(i) == ' 28Si' ) THEN
        i_28Si            = i
      END IF 
      IF (      a_name(i) == ' 32S ' ) THEN
        i_32S             = i
      END IF 
      IF (      a_name(i) == ' 36Ar' ) THEN
        i_36Ar            = i
      END IF 
      IF (      a_name(i) == ' 40Ca' ) THEN
        i_40Ca            = i
      END IF 
      IF (      a_name(i) == ' 44Ti' ) THEN
        i_44Ti            = i
      END IF 
      IF (      a_name(i) == ' 48Cr' ) THEN
        i_48Cr            = i
      END IF 
      IF (      a_name(i) == ' 52Fe' ) THEN
        i_52Fe            = i
      END IF 
      IF (      a_name(i) == ' 56Ni' ) THEN
        i_56Ni            = i
      END IF 
      IF (      a_name(i) == ' 60Zn' ) THEN
        i_60Zn            = i
      END IF 
      IF (      a_name(i) == '  n  ' ) THEN
        i_n               = i
      END IF 
      IF (      a_name(i) == '  p  ' ) THEN
        i_p               = i
      END IF 
      IF (      a_name(i) == '  4He' ) THEN
        i_4He             = i
      END IF 
    END DO ! i = 1,nuc_number
  END IF ! first_c

!-----------------------------------------------------------------------
!  NSE-nonNSE boundary
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    DO k = kmin,kmax
      DO i = imin,imax
        IF ( nseMD(i,j,k) == 0 ) THEN
          i_nse(j,k)      = i
          r_nse(j,k)      = x_e_in(i)
          EXIT
        END IF ! nseMD(i,j,k) == 0
      END DO ! i = imin,imax
    END DO ! k = kmin,kmax
  END DO ! j = jmin,jmax

!-----------------------------------------------------------------------
!  0.1 mass 16O fraction boundary
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    DO k = kmin,kmax
      DO i = i_nse(j,k),imax
        IF ( xnMD(i,i_16O,j,k) > 0.1d0 ) THEN
          i_O1u           = MAX( i    , 2 )
          i_O1l           = MAX( i - 1, 1 )
          EXIT
        END IF ! xnMD(i,i_16O,j) > 0.1d0
      END DO ! i = i_nse(j,k),imax
      r_O1(j,k)           = rinterp( x_e_in(i_O1u), x_e_in(i_O1l), xnMD(i_O1u,i_16O,j,k), 0.1d0, xnMD(i_O1l,i_16O,j,k) )
    END DO ! k = kmin,kmax
  END DO ! j = jmin,jma

!-----------------------------------------------------------------------
!  0.5 mass 16O fraction boundary
!-----------------------------------------------------------------------

y: DO j = jmin,jmax
z:   DO k = kmin,kmax
      DO i = i_nse(j,k),imax
        IF ( xnMD(i,i_16O,j,k) > x_O ) THEN
          i_xOu           = MAX( i    , 2 )
          i_xOl           = MAX( i - 1, 1 )
          EXIT
        END IF ! xnMD(i,i_16O,j) > x_O
        IF ( i == imax ) THEN
          r_xO(:,:)       = zero
          EXIT z
        END IF ! i == imax
      END DO ! i = i_nse(j,k),imax
      r_xO(j,k)           = rinterp( x_e_in(i_xOu), x_e_in(i_xOl), xnMD(i_xOu,i_16O,j,k), x_O, xnMD(i_xOl,i_16O,j,k) )
    END DO z
  END DO y

!-----------------------------------------------------------------------
!  Dominant species
!
!    i_comp = 1 : free neutrons and protons
!    i_comp = 2 : helium
!    i_comp = 3 : oxygen, neon and magnesium
!    i_comp = 4 : silicon, sulfur and argon
!    i_comp = 5 : Ni-like elements
!    i_comp = 6 : auxiliary heavy nucleus
!    i_comp = 7 : oxygen, neon and magnesium exceeding 10 percent, but not dominant
!-----------------------------------------------------------------------

  r_comp                  = zero
  DO i = imin,imax
    DO j = jmin,jmax
      DO k = kmin,kmax
        IF ( nseMD(i,j,k) == 0 ) THEN
          r_comp(i,1,j,k) = xnMD(i,i_n,j,k) + xnMD(i,i_p,j,k)
          r_comp(i,2,j,k) = xnMD(i,i_4He ,j,k)
          r_comp(i,3,j,k) = xnMD(i,i_16O ,j,k) + xnMD(i,i_20Ne,j,k) + xnMD(i,i_24Mg,j,k)
          r_comp(i,4,j,k) = xnMD(i,i_28Si,j,k) + xnMD(i,i_32S ,j,k) + xnMD(i,i_36Ar,j,k)
          r_comp(i,5,j,k) = xnMD(i,i_40Ca,j,k) + xnMD(i,i_44Ti,j,k) + xnMD(i,i_48Cr,j,k) &
&                         + xnMD(i,i_52Fe,j,k) + xnMD(i,i_56Ni,j,k) + xnMD(i,i_60Zn,j,k)
          r_comp(i,6,j,k) = xnMD(i,nuc_number+1,j,k)
          i_comp_mx       = MAXLOC(r_comp(i,:,j,k))
          i_comp(i,j,k)   = i_comp_mx(1)
          IF ( i_comp(i,j,k) /= 3  .and.  r_comp(i,3,j,k) > 0.1d0 ) i_comp(i,j,k) = 7
        ELSE ! nseMD(i,j,k) /= 0
          r_comp(i,1,j,k) = aesvMD(i+1,7,j,k) + aesvMD(i+1,8,j,k)
          r_comp(i,6,j,k) = aesvMD(i+1,9,j,k)
          r_comp(i,2,j,k) = one - r_comp(i,1,j,k) - r_comp(i,6,j,k)
          i_comp_mx       = MAXLOC(r_comp(i,:,j,k))
          i_comp(i,j,k)   = i_comp_mx(1)
        END IF ! nseMD(i,j,k) == 0
      END DO ! j = jmin,jmax
    END DO ! k = kmin,kmax
  END DO ! i = imin,imax

!-----------------------------------------------------------------------
!
!               \\\\\ TRANSFER VARIABLES TO EDIT /////
!
!-----------------------------------------------------------------------

CALL edit_MD_exec( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, nnu,     &
& time, t_tb, x_c_in, y_in, z_in, rhoMD, uMD, vMD, wMD, pMD, eMD, sMD,      &
& aMD, yeMD, v_csoundMD, e_nu_MD, f_nu_MD, lumMD, e_rmsMD, r_nse, r_O1,     &
& r_xO, r_shock, r_shock_mn, r_shock_mx, rsphere_mean, dsphere_mean,        &
& tsphere_mean, msphere_mean, esphere_mean, r_gain, tau_advMD,              &
& tau_heat_nuMD, tau_heat_nucMD, n_growMD, wBVMD, ylMD, i_comp, dudt_nucMD, &
& dudt_nuMD, grav_x_cMD, grav_y_cMD, grav_z_cMD, l_editMDu, l_editMDv,      &
& l_editMDw, l_editMDs, l_editMDd, l_editMDe, l_editMDp, l_editMDenu,       &
& l_editMDfnu, l_editMDa, l_editMDx, l_editMDye, l_editMDcm, l_editMDnu,    &
& l_editMDnc, l_editMDnl, l_editMDne, l_editMDgx, l_editMDgy, l_editMDgz,   &
& l_editMDBVw, l_editMDyl, l_MDedit )

!-----------------------------------------------------------------------
!               ||||| MYID == 0 ONLY ENDS HERE |||||
!-----------------------------------------------------------------------

END IF ! myid == 0

CALL MPI_Barrier( MPI_COMM_WORLD, ierr )

RETURN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CONTAINS
  REAL (KIND=double) FUNCTION mean(q1,q2,x1,x2)

  REAL (KIND=double) :: q1
  REAL (KIND=double) :: q2
  REAL (KIND=double) :: x1
  REAL (KIND=double) :: x2

  mean                 = ( q1 * x1 + q2 * x2 )/( q1 + q2 )

END FUNCTION mean

  REAL (KIND=double) FUNCTION rinterp(a,b,x,y,z)

  REAL (KIND=double) :: a
  REAL (KIND=double) :: b
  REAL (KIND=double) :: x
  REAL (KIND=double) :: y
  REAL (KIND=double) :: z

  rinterp      = b + ( a - b ) * ( y - z )/( x - z )

END FUNCTION rinterp

END SUBROUTINE edit_MD_inout

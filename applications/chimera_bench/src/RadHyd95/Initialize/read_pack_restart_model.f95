SUBROUTINE read_pack_restart_model( nreadp, nprint, iskip, nx, nez, nnu,   &
& nnc, ij_ray_dim, ik_ray_dim, i_model_data, d_model_data1, d_model_data2, &
& d_model_data3, i_nuc_data, d_nuc_data, d_psi_data1, d_psi_data2,         &
& d_psi_data3, d_psi_data4, d_psi_data5, nuc_number, nrst )
!-----------------------------------------------------------------------
!
!    File:         read_pack_restart_model
!    Module:       read_pack_restart_model
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the restart model configuration and nuclear abundances
!       and to pack the data into an integer and a real*8 array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp       : unit number from which to read
!  nprint       : unit number from which to print
!  iskip        : echo data read flag
!  nx           : radial array dimension
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  nnc          : nuclear specie array dimension
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  nuc_number   : number of nuclear species (not counting representative heavy nucleus)
!  nrst         : cycle number at start or restart
!
!    Output arguments:
!  i_model_data  : integer array of initial model data
!  d_model_data1 : 64 bit real array of initial model data
!  d_model_data2 : 64 bit real array of initial model data
!  d_model_data3 : 64 bit real array of initial model data
!  i_nuc_data    : integer array of nuclear keys
!  d_nuc_data    : real*8 array of nuclear keys
!  d_psi_data1   : 64 bit real array of neutrino data
!  d_psi_data2   : 64 bit real array of neutrino data
!  d_psi_data3   : 64 bit real array of neutrino data
!  d_psi_data4   : 64 bit real array of neutrino data
!  d_psi_data5   : 64 bit real array of neutrino data
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc
USE numerical_module, ONLY : zero, epsilon

USE edit_module, ONLY : nlog

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc           ! nuclear specie array dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nuc_number    ! number of nuclear species (not counting representative heavy nucleus)
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(2)                                       :: i_model_data ! integer array of initial model data

REAL(KIND=double), INTENT(out), DIMENSION(7,nx,ij_ray_dim,ik_ray_dim)                :: d_model_data1 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx)                                      :: d_model_data2 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx+1)                                    :: d_model_data3 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(10+4*nnc+(nnc+4)*nx,ij_ray_dim,ik_ray_dim) :: d_nuc_data    ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(8,nez,nnu,ij_ray_dim,ik_ray_dim)           :: d_psi_data1   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(2,nx,nez,nnu,ij_ray_dim,ik_ray_dim)        :: d_psi_data2   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(2,nx,nnu,ij_ray_dim,ik_ray_dim)            :: d_psi_data3   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(2,nnu,ij_ray_dim,ik_ray_dim)               :: d_psi_data4   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(out), DIMENSION(2,ij_ray_dim,ik_ray_dim)                   :: d_psi_data5   ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(inout), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim)                       :: i_nuc_data   ! integer array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: i             ! do index
INTEGER                          :: j             ! radial zone index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status
INTEGER                          :: ij_ray        ! radial ray j index
INTEGER                          :: ik_ray        ! radial ray k index
INTEGER                          :: iadjst        ! mass fraction normalization flag

REAL(KIND=double)                :: xn_tot        ! sum of themass fractions in a zone

!-----------------------------------------------------------------------
!  Radial mass zoning
!-----------------------------------------------------------------------
!  imin   : the innermost unshifted radial zone
!
!  imax   : the outermost unshifted radial zone
!
!  jr_min : the innermost MGFLD shifted radial zone
!
!  jr_max : the outermost MGFLD shifted radial zone
!-----------------------------------------------------------------------

INTEGER                                            :: imin
INTEGER                                            :: imax
INTEGER                                            :: jr_min
INTEGER                                            :: jr_max

!-----------------------------------------------------------------------
!  Thermodynamic state
!-----------------------------------------------------------------------
!  rho(j,ij_ray,ik_ray) : the density of radial zone j at timestep m (current
!   time) (g/cm**3).
!
!  t(j,ij_ray,ik_ray)   : the temperature of radial zone j at timestep m (K).
!
!  ye(j,ij_ray,ik_ray)  : the electron fraction of radial zone j at timestep m.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: t
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: ye

!-----------------------------------------------------------------------
!  Mechanical state
!-----------------------------------------------------------------------
!  u(j,ij_ray,ik_ray)  : the zone averaged x-velocity of radial zone j at
!   timestep m (cm/s).
!
!  v(j,ij_ray,ik_ray)  : the zone averaged y-velocity of radial zone j at
!   timestep m (cm/s).
!
!  e(j,ij_ray,ik_ray)  : the zone averaged z-velocity of radial zone j at
!   timestep m (cm/s).
!
!  r(j)        : the radius of outer boundary of radial zone j at
!   timestep m (cm).
!
!  dr(j)       : r(j) - r(j-1).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: u
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: v
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: w
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dr

!-----------------------------------------------------------------------
!  Radiation Variables
!-----------------------------------------------------------------------
!  psi0(j,k,n,ij_ray,ik_ray) : the zero moment of the occupation distribution
!   for neutrinos at the midpoint of radial zone j, of energy zone k,
!   and of type n.
!
!  e_nu_c_bar(j)       : angular averaged neutrino energy density
!   (ergs cm^{-3}).
!
!  f_nu_e_bar(j)       : angular averaged neutrino energy flux
!   (ergs cm^{-2} s^{-1}).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: psi0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: e_nu_c_bar
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: f_nu_e_bar

!-----------------------------------------------------------------------
!  Transport parameters derived from psi0 and psi1
!-----------------------------------------------------------------------
!  dnurad(j,k,n,ij_ray,ik_ray) : the total number of n-neutrinos of
!   energy zone k per unit energy that have crossed the outer boundary
!   radial zone j (/MeV).
!
!  unukrad(k,n,ij_ray,ik_ray)  : the cumulative energy emitted from the
!   core by n-type neutrinos of energy zone k (ergs).
!
!  unujrad(j,n,ij_ray,ik_ray)  : the cumulative energy transported across
!   the outer boundary of radial zone j by n-type neutrinos (ergs).
!
!  nnukrad(k,n,ij_ray,ik_ray)  : the cumulative number of n-type
!   neutrinos of energy k emitted by the core.
!
!  nnujrad(j,n,ij_ray,ik_ray)  : the net number of n-type neutrinos
!   transported across the outer boundary of radial zone j.
!
!  e_rad(ij_ray,ik_ray)        : the cumulative material energy entering
!   (negative) or leaving (positive) the grid (ergs).
!
!  unurad(n,ij_ray,ik_ray)     : the cumulative energy emitted from the
!   core in n-type neutrinos (ergs).
!
!  elec_rad(ij_ray,ik_ray)     : the net number of electrons advected in
!   (negative) or out (positive) of the grid
!
!  nnurad(n,ij_ray,ik_ray)     : the cumulative number of n-type
!   neutrinos emitted by the core.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)  :: dnurad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unukrad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: unujrad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: nnukrad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)    :: nnujrad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: e_rad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: unurad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)        :: elec_rad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)      :: nnurad

!-----------------------------------------------------------------------
!  NSE flag
!-----------------------------------------------------------------------
!  nse : nuclear statistical equilibrium flag
!
!     nse(i,ij_ray,ik_ray) = 0 : material not in nuclear statistical
!      equilibrium; nuclear reaction network must be turned on to evolve
!      the matter composition.
!     nse(i,ij_ray,ik_ray) = 1 : material in nuclear statistical
!      equilibrium; nuclear reaction network turned off.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nse

!-----------------------------------------------------------------------
!  Abundance parameters
!-----------------------------------------------------------------------
!  xn(j,i,ij_ray,ik_ray)       : mass fraction of the ith nucleus.
!
!  uburn(j,ij_ray,ik_ray)      : cumulative energy generated in zone j by nuclear
!   reactions (ergs/gm).
!
!  be_nuc_rep(j,ij_ray,ik_ray) : binding energy of the representative heavy
!   nucleus (MeV).
!
!  a_nuc_rep(j,ij_ray,ik_ray)  : mass number of the representative heavy nucleus.
!
!  z_nuc_rep(j,ij_ray,ik_ray)  : charge number of the representative heavy
!   nucleus.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: xn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: uburn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: be_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: a_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: z_nuc_rep

!-----------------------------------------------------------------------
!  Parameters for nurad plot edits
!-----------------------------------------------------------------------
!  nu_r(k,n,ij_ray,ik_ray)     : the number of neutrinos of energy group
!   k radiated across r_nurad in time dtnuradplot.
!
!  nu_rt(k,n,ij_ray,ik_ray)    : the cumulative number of neutrinos of
!   energy group k radiated across r_nurad.
!
!  nu_rho(k,n,ij_ray,ik_ray)   : the number of neutrinos of energy group
!   k radiated  across rho_nurad in time dtnuradplot.
!
!  nu_rhot(k,n,ij_ray,ik_ray)  : the cumulative number of neutrinos of
!   energy group k radiated across rho_nurad.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: nu_r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: nu_rt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: nu_rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: nu_rhot

!-----------------------------------------------------------------------
!  Parameters for neutrino distribution plot edits
!-----------------------------------------------------------------------
!  psi0dat(k,n,ij_ray,ik_ray) : the time integrated psi0, to be time
!   averaged on dumping.
!
!  psi1dat(k,n,ij_ray,ik_ray) : the time integrated psi1, to be time
!   averaged on dumping.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: psi0dat
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)  :: psi1dat

!-----------------------------------------------------------------------
!  Energy offset
!-----------------------------------------------------------------------
!  duesrc(j,ij_ray,ik_ray) : cumulative energy offsets generated by
!   regridding the EOS [ergs g^{-1}].
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)   :: duesrc

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in read_pack_initial_model')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_pack_initial_model')
 3001 FORMAT (' jr_max + 1 =',i4,' > nx =',i4)
 3003 FORMAT (' jr_max,n_proc =',2i6,' MOD( jr_max - 1 ,n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 )')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (rho(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi0(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_nu_c_bar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_nu_c_bar'; WRITE  (nlog,1001) var_name; END IF
ALLOCATE (f_nu_e_bar(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f_nu_e_bar'; WRITE  (nlog,1001) var_name; END IF

ALLOCATE (dnurad(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dnurad    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unukrad(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unukrad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unujrad(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujrad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnukrad(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnukrad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnujrad(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnujrad   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e_rad(ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_rad     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unurad(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unurad    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (elec_rad(ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'elec_rad  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nnurad(nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnurad    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (nse(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (xn(nx,nnc,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (be_nuc_rep(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (a_nuc_rep(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_nuc_rep(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uburn(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uburn     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (nu_r(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_r      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_rt(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rt     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_rho(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rho    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nu_rhot(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rhot   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi0dat(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0dat   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1dat(nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1dat   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (duesrc(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'duesrc    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

imin                   = 0
imax                   = 0
jr_min                 = 0
jr_max                 = 0
nse                    = 0

rho                    = zero
t                      = zero
ye                     = zero
u                      = zero
r                      = zero
dr                     = zero
psi0                   = zero

dnurad                 = zero
unukrad                = zero
unujrad                = zero
e_rad                  = zero
unurad                 = zero
nnukrad                = zero
nnujrad                = zero
elec_rad               = zero
nnurad                 = zero

xn                     = zero
be_nuc_rep             = zero
a_nuc_rep              = zero
z_nuc_rep              = zero
uburn                  = zero

!-----------------------------------------------------------------------
!
!                \\\\\ READ INITIAL MODEL DATA /////
!
!-----------------------------------------------------------------------

REWIND (nreadp)

!-----------------------------------------------------------------------
!  Radial index bounds for unshifted radial arrays
!-----------------------------------------------------------------------

READ (nreadp) imin
READ (nreadp) imax

jr_min                 = imin + 1
jr_max                 = imax + 1

!-----------------------------------------------------------------------
!  Independent thermodynamic variables
!-----------------------------------------------------------------------

READ (nreadp) rho
READ (nreadp) t
READ (nreadp) ye

!-----------------------------------------------------------------------
!  Independent mechanical variables
!-----------------------------------------------------------------------

READ (nreadp) u
READ (nreadp) v
READ (nreadp) w
READ (nreadp) dr
READ (nreadp) r

!-----------------------------------------------------------------------
!  Independent radiation variables and bookkeeping arrays
!-----------------------------------------------------------------------

READ (nreadp) psi0

!-----------------------------------------------------------------------
!  unucr
!-----------------------------------------------------------------------

READ (nreadp) dnurad
READ (nreadp) unukrad
READ (nreadp) unujrad
READ (nreadp) e_rad
READ (nreadp) unurad

!-----------------------------------------------------------------------
!  nnucr
!-----------------------------------------------------------------------

READ (nreadp) nnukrad
READ (nreadp) nnujrad
READ (nreadp) elec_rad
READ (nreadp) nnurad

!-----------------------------------------------------------------------
!  Angular averaged neutrino energy and fluxes
!-----------------------------------------------------------------------

READ (nreadp) e_nu_c_bar
READ (nreadp) f_nu_e_bar

!-----------------------------------------------------------------------
!  Net number of neutrinos radiated from density rho_nurad and radius
!   r_nurad
!-----------------------------------------------------------------------

READ (nreadp) nu_r
READ (nreadp) nu_rt
READ (nreadp) nu_rho
READ (nreadp) nu_rhot

!-----------------------------------------------------------------------
!  Time integrated psi0 and psi1
!-----------------------------------------------------------------------

READ (nreadp) psi0dat
READ (nreadp) psi1dat

!-----------------------------------------------------------------------
!  nse - non-bse flag
!-----------------------------------------------------------------------

READ (nreadp) nse

!-----------------------------------------------------------------------
!  Nuclear abundances
!-----------------------------------------------------------------------

READ (nreadp) xn

!-----------------------------------------------------------------------
!  Auxiliary heavy nucleus
!-----------------------------------------------------------------------

READ (nreadp) a_nuc_rep
READ (nreadp) z_nuc_rep
READ (nreadp) be_nuc_rep

!-----------------------------------------------------------------------
!  Nuclear energy released
!-----------------------------------------------------------------------

READ (nreadp) uburn

!-----------------------------------------------------------------------
!  Energy offsets
!-----------------------------------------------------------------------

READ (nreadp, IOSTAT = istat) duesrc

!-----------------------------------------------------------------------
!
!             \\\\\ CHECK DATA EXTENTS FOR CONSISTENCY /////
!
!  jr_max must be compatible with n_proc, so the y_arrays and z-arrays
!   can be distributed uniformly over the processors.
!
!-----------------------------------------------------------------------

IF ( jr_max + 1 > nx ) THEN
  WRITE (nlog,3001) jr_max+1,nx
  WRITE (nprint,3001) jr_max,nx
  STOP
END IF ! jr_max + 1 > nx

IF ( MOD( jr_max - 1 , n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 ) /= 0 ) THEN
  WRITE (nlog,3003) jr_max,n_proc
  WRITE (nprint,3003) jr_max,n_proc
  STOP
END IF ! MOD( jr_max - 1 ,n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 )

!-----------------------------------------------------------------------
!
!                \\\\\ ADJUST NUCLEAR ABUNDANCES /////
!
!-----------------------------------------------------------------------

iadjst                         = i_nuc_data(4,1,1)

IF ( iadjst == 1 ) THEN
  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim
      DO j = 1,nx
        IF ( nse(j,ij_ray,ik_ray) == 0 ) THEN
          xn_tot               = zero
          DO n = 1,nuc_number+1
            xn_tot             = xn_tot + xn(j,n,ij_ray,ik_ray)
          END DO
          DO n = 1,nuc_number+1
            xn(j,n,ij_ray,ik_ray) = xn(j,n,ij_ray,ik_ray)/( xn_tot + epsilon )
          END DO
        END IF ! nse(j,ij_ray,ik_ray) == 0
      END DO ! j
    END DO ! ij_ray
  END DO ! ik_ray
END IF ! iadjst == 1

!-----------------------------------------------------------------------
!
!                \\\\\ PACK INITIAL MODEL DATA /////
!
!-----------------------------------------------------------------------

i_model_data(1)                = imin
i_model_data(2)                = imax

d_model_data1(1,:,:,:)         = rho   (:,:,:)
d_model_data1(2,:,:,:)         = t     (:,:,:)
d_model_data1(3,:,:,:)         = ye    (:,:,:)
d_model_data1(4,:,:,:)         = u     (:,:,:)
d_model_data1(5,:,:,:)         = v     (:,:,:)
d_model_data1(6,:,:,:)         = w     (:,:,:)
d_model_data1(7,:,:,:)         = duesrc(:,:,:)

d_model_data2(1,:)             = dr        (:)
d_model_data2(2,:)             = e_nu_c_bar(:)

d_model_data3(1,:)             = r         (:)
d_model_data3(2,:)             = f_nu_e_bar(:)

d_psi_data1(1,:,:,:,:)         = nu_r   (:,:,:,:)
d_psi_data1(2,:,:,:,:)         = nu_rt  (:,:,:,:)
d_psi_data1(3,:,:,:,:)         = nu_rho (:,:,:,:)
d_psi_data1(4,:,:,:,:)         = nu_rhot(:,:,:,:)
d_psi_data1(5,:,:,:,:)         = psi0dat(:,:,:,:)
d_psi_data1(6,:,:,:,:)         = psi1dat(:,:,:,:)
d_psi_data1(7,:,:,:,:)         = unukrad(:,:,:,:)
d_psi_data1(8,:,:,:,:)         = nnukrad(:,:,:,:)

d_psi_data2(1,:,:,:,:,:)       = psi0  (:,:,:,:,:)
d_psi_data2(2,:,:,:,:,:)       = dnurad(:,:,:,:,:)

d_psi_data3(1,:,:,:,:)         = unujrad(:,:,:,:)
d_psi_data3(2,:,:,:,:)         = nnujrad(:,:,:,:)

d_psi_data4(1,:,:,:)           = unurad(:,:,:)
d_psi_data4(2,:,:,:)           = nnurad(:,:,:)

d_psi_data5(1,:,:)             = e_rad   (:,:)
d_psi_data5(2,:,:)             = elec_rad(:,:)

!-----------------------------------------------------------------------
!
!                \\\\\ PACK NUCLEAR ABUNDANCE DATA /////
!
!-----------------------------------------------------------------------

i_nuc_data(11:nx+10,:,:)       = nse(1:nx,:,:)

DO n = 1,nnc
  DO i = 1,nx
    d_nuc_data(10+4*nnc+(n-1)*nx+i,:,:) = xn(i,n,:,:)
  END DO ! i
END DO ! n

DO i = 1,nx
  d_nuc_data(10+4*nnc+(nnc+0)*nx+i,:,:) = a_nuc_rep(i,:,:)
  d_nuc_data(10+4*nnc+(nnc+1)*nx+i,:,:) = z_nuc_rep(i,:,:)
  d_nuc_data(10+4*nnc+(nnc+2)*nx+i,:,:) = be_nuc_rep(i,:,:)
  d_nuc_data(10+4*nnc+(nnc+3)*nx+i,:,:) = uburn(i,:,:)
END DO ! i

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (rho, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (u, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (v, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (w, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (r, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (psi0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (e_nu_c_bar, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_nu_c_bar'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (f_nu_e_bar, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'f_nu_e_bar'; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dnurad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dnurad    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (unukrad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unukrad   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (unujrad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unujrad   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (nnukrad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnukrad   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (nnujrad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnujrad   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (e_rad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_rad     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (unurad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unurad    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (elec_rad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'elec_rad  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (nnurad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnurad    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (nse, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nse       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (xn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (be_nuc_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'be_nuc_rep'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (a_nuc_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'a_nuc_rep '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (z_nuc_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_nuc_rep '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (uburn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uburn     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (nu_r, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_r      '; WRITE (nlog,1001) var_name; END IF
DEALLOCATE (nu_rt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rt     '; WRITE (nlog,1001) var_name; END IF
DEALLOCATE (nu_rho, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rho    '; WRITE (nlog,1001) var_name; END IF
DEALLOCATE (nu_rhot, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nu_rhot   '; WRITE (nlog,1001) var_name; END IF

DEALLOCATE (psi0dat, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0dat   '; WRITE (nlog,1001) var_name; END IF
DEALLOCATE (psi1dat, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1dat   '; WRITE (nlog,1001) var_name; END IF

DEALLOCATE (duesrc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'duesrc    '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE read_pack_restart_model

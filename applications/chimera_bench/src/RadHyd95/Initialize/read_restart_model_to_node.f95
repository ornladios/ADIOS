SUBROUTINE read_restart_model_to_node( nreadp, nprint, iskip, nx, nez, &
& nnu, nnc, nrst, ij_ray_dim, ik_ray_dim, i_nuc_data, nuc_number,      &
& i_model_data )
!-----------------------------------------------------------------------
!
!    File:         read_restart_model_to_node
!    Module:       read_restart_model_to_node
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the initial model directly to a specified processor.
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
!  nrst         : cycle number at start or restart
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  nuc_number   : number of nuclear species (not counting representative heavy nucleus)
!
!    Output arguments:
!  i_model_data : minimum and maximum of active radial ray zone centers
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module, nucbrn_module,
!  nu_dist_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc
USE numerical_module, ONLY : zero, half, epsilon

USE edit_module, ONLY : nlog, nu_r, nu_rt, nu_rho, nu_rhot, psi0dat, psi1dat
USE eos_snc_x_module, ONLY : nse_e=>nse, duesrc
USE mdl_cnfg_module, ONLY : jr_min, jr_max
USE nucbrn_module, ONLY : nse_n=>nse
USE nu_dist_module, ONLY : dnurad, unukrad, unujrad, nnukrad, nnujrad, &
& e_rad, unurad, elec_rad, nnurad
USE radial_ray_module, ONLY : rho_c, t_c, ye_c, u_c, v_c, w_c, dx_ci, &
& x_ei, x_ci, psi0_c, xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, &
& uburn_c, nse_c, e_nu_c_bar, f_nu_e_bar

IMPLICIT none
SAVE

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
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nuc_number    ! number of nuclear species (not counting representative heavy nucleus)

INTEGER, INTENT(in), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim) :: i_nuc_data ! integer array of edit keys

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(2) :: i_model_data  ! integer array of initial model data

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

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
!-----------------------------------------------------------------------

INTEGER                                              :: imin
INTEGER                                              :: imax

!-----------------------------------------------------------------------
!  Thermodynamic state
!-----------------------------------------------------------------------
!  rho(j) : the density of radial zone j at timestep m (current time)
!   (g/cm**3).
!
!  t(j)   : the temperature of radial zone j at timestep m (K).
!
!  ye(j)  : the electron fraction of radial zone j at timestep m.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: t
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: ye

!-----------------------------------------------------------------------
!  Mechanical state
!-----------------------------------------------------------------------
!  u(j)  : the zone averaged x-velocity of radial zone j at timestep m
!   (cm/s).
!
!  v(j)  : the zone averaged y-velocity of radial zone j at timestep m
!   (cm/s).
!
!  e(j)  : the zone averaged z-velocity of radial zone j at timestep m
!   (cm/s).
!
!  r(j)  : the radius of outer boundary of radial zone j at timestep m
!   (cm).
!
!  dr(j) : r(j) - r(j-1).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: u
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: v
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: w
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dr

!-----------------------------------------------------------------------
!  Radiation Variables
!-----------------------------------------------------------------------
!  psi0(j,k,n) : the zero moment of the occupation distribution for
!   neutrinos at the midpoint of radial zone j, of energy zone k, and of
!   type n.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: psi0

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
!  uburn(j,ij_ray,ik_ray)      : cumulative energy generated in zone j
!   by nuclear reactions (ergs/gm).
!
!  be_nuc_rep(j,ij_ray,ik_ray) : binding energy of the representative
!   heavy nucleus (MeV).
!
!  a_nuc_rep(j,ij_ray,ik_ray)  : mass number of the representative heavy
!   nucleus.
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
!  Energy offset
!-----------------------------------------------------------------------
!  duesrc_r(j,i_ray) : cumulative energy offsets generated by regridding
!   the EOS [ergs g^{-1}].
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)   :: duesrc_r

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in read_pack_initial_model')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_pack_initial_model')
 3001 FORMAT (' jr_max + 1 =',i4,' > nx =',i4)
 3003 FORMAT (' jr_max,n_proc =',2i6,' MOD( jr_max - 1 ,n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 ) /= 0')

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
ALLOCATE (duesrc_r(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'duesrc_r  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!        Initialize
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

i_model_data(1)        = imin
i_model_data(2)        = imax

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

READ (nreadp, IOSTAT = istat) duesrc_r

!-----------------------------------------------------------------------
!
!             \\\\\ CHECK DATA EXTENTS FOR CONSISTENCY /////
!
!        jr_max must be compatible with n_proc, so the y_arrays and
!         z-arrays can be distributed uniformly over the processors.
!
!-----------------------------------------------------------------------

IF ( jr_max + 1 > nx ) THEN
  WRITE (nlog,3001) jr_max+1,nx
  WRITE (nprint,3001) jr_max+1,nx
  STOP
END IF ! jr_max + 1 > nx

IF ( MOD( jr_max - 1 ,n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 ) /= 0 ) THEN
  WRITE (nlog,3003) jr_max,n_proc
  WRITE (nprint,3003) jr_max,n_proc
  STOP
END IF ! MOD( jr_max - 1 ,n_proc ) /= 0  .and.  MOD( n_proc, jr_max - 1 )

!-----------------------------------------------------------------------
!
!                \\\\\ ADJUST NUCLEAR ABUNDANCES /////
!
!-----------------------------------------------------------------------

iadjst                 = i_nuc_data(4,1,1)

IF ( iadjst == 1 ) THEN
  DO ik_ray = 1,ik_ray_dim
    DO ij_ray = 1,ij_ray_dim
      DO j = 1,nx
        IF ( nse(j,ij_ray,ik_ray) == 0 ) THEN
          xn_tot         = zero
          DO n = 1,nuc_number+1
            xn_tot       = xn_tot + xn(j,n,ij_ray,ik_ray)
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
!        \\\\\ TRANSFER INITIAL DATA TO RADIAL_RAY_MODULE /////
!
!-----------------------------------------------------------------------


rho_c                   = rho
t_c                     = t
ye_c                    = ye
u_c                     = u
v_c                     = v
w_c                     = w
x_ei                    = r
dx_ci                   = dr
x_ci(1:nx)              = x_ei(1:nx) + half * dx_ci(1:nx)
psi0_c                  = psi0

nse_c                   = nse
nse_e(2:nx,:,:)         = nse(1:nx-1,:,:)
nse_n(2:nx)             = nse(1:nx-1,1,1)

xn_c                    = xn
a_nuc_rep_c             = a_nuc_rep
z_nuc_rep_c             = z_nuc_rep
be_nuc_rep_c            = be_nuc_rep
uburn_c                 = uburn
duesrc                  = duesrc_r

!-----------------------------------------------------------------------
!        Deallocate arrays
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
DEALLOCATE (duesrc_r, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'duesrc_r  '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE read_restart_model_to_node

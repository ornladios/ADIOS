SUBROUTINE setup_restart_hydro( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, nx, ny, nz, nnc, rhop, tp, yep, ye_ip, rp, drp, up, xnp, &
& be_nuc_repp, a_nuc_repp, z_nuc_repp, aesv_c )
!-----------------------------------------------------------------------
!
!    File:         setup_restart_hydro
!    Module:       setup_restart_hydro
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!    Purpose:
!      To compute, from a restart file of problem configuration,
!       the quantities needed to continue executing the simulation.
!       In particular,
!
!    (1) the equation of state tables are initialized;
!    (2) radii, etc are computed;
!    (3) variable arrays at time n-1 and time n+1 are filled;
!    (4) problem modifications are made (through subroutine 'pblmst1')
!
!    Subprograms called:
!  pblmst1       : option to change problem before neutrino quantities are calculated
!  esrgnz_x      : loads EOS tables
!  eqstz_x       : interpolates EOS quantities
!  gammaz_x      : computes the EOS gammas
!
!    Input arguments:
!  imin          : inner x-zone array index
!  imax          : outer x-zone array index
!  ij_ray        : j-index of a radial ray
!  ik_ray        : k-index of a radial ray
!  ij_ray_dim    : number of y-zones on a processor before swapping
!  ik_ray_dim    : number of z-zones on a processor before swapping
!  nx            : x_array extent
!  ny            : y_array extent
!  nz            : z_array extent
!  nnc           : neutrino abundance array extent
!  rhop          : density (cm^{-3})
!  tp            : temperature (MeV)
!  yep           : electron fraction
!  ye_ip         : initial electron fraction
!  rp            : radial zone radii (cm)
!  up            : radial velocity of zone (cm)
!  xnp           : initial mass fractions
!  be_nuc_repp   : binding energy of mean heavy nucleus
!  a_nuc_repp    : mass number of mean heavy nucleus
!  z_nuc_repp    : charge number of mean heavy nucleus
!
!    Output arguments:
!  aesv_c        : intwrpolated EOS quantities
!  yep           : electron fraction
!  ye_ip         : initial electron fraction
!  up            : radial velocity of zone (cm)
!
!    Include files:
!  kind_module, numerical_module
!  boundary_module, convect_module, cycle_module, edit_module,
!  eos_snc_x_module, mdl_cnfg_module, nucbrn_module, shock_module
!     
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, epsilon

USE boundary_module, ONLY : iubcjmn, iubcjmx, ubcjmn, ubcjmx, iuset, uset, &
& ipbnd, pbound
USE convect_module
USE cycle_module, ONLY : ncycle, nrst
USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY : aesv, xn, be_nuc_rep, a_nuc_rep, z_nuc_rep, nse
USE mdl_cnfg_module, ONLY : jr_min, jr_max, jmr, r, rhor, rho, tr, t, yer, &
& ye, ye0, u, dr, rho0, p0
USE nucbrn_module, ONLY : a_nuc, z_nuc, nuc_number, a_name
USE shock_module, ONLY : q0_x, q0_y, q0_z

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin            ! inner x_array index
INTEGER, INTENT(in)              :: imax            ! outer x_array index

INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim      ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of z-zones on a processor before swapping

INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: ny              ! y-array extent
INTEGER, INTENT(in)              :: nz              ! z-array extent
INTEGER, INTENT(in)              :: nnc             ! composition array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: rhop        ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: tp          ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                         :: rp          ! radial zone radii (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)                           :: drp         ! radial zone thickness (cm)

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: be_nuc_repp ! binding energy of auxiliary nucleus (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: a_nuc_repp  ! mass number of auxiliary nucleus
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: z_nuc_repp  ! charge number of auxiliary nucleus

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,12,ij_ray_dim,ik_ray_dim) :: aesv_c      ! unsifted zone-edged lapse function
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: ye_ip       ! electron fraction

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: up       ! radial velocity of zone (cm s^{-1})
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: yep      ! electron fraction
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim) :: xnp      ! mass fractions

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first_c = .true.

INTEGER                          :: i               ! composition index
INTEGER                          :: j               ! radial zone index
INTEGER                          :: jr_maxp         ! jr_max + 1

INTEGER                          :: i_n             ! neutron abundance index
INTEGER                          :: i_p             ! proton abundance index
INTEGER                          :: i_4He           ! 4He abundance index

REAL(KIND=double)                :: z_tot           ! total proton number (used for computing ye)
REAL(KIND=double)                :: a_tot           ! total mass number (used for computing ye)

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Variable initialization complete in setup_restart_hydro')
  103 FORMAT (' Velocity boundary conditions set (if option selected) in setup_restart_hydro')
  105 FORMAT (' Variable transfer complete in setup_restart_hydro')
  107 FORMAT (' pblmst1 called in setup_restart_hydro')
  109 FORMAT (' Equation of state tables are being filled in setup_restart_hydro')
  111 FORMAT (' Equation of state quantities are being interpolated in setup_restart_hydro')
  113 FORMAT (' Adiabatic exponents are being computed in setup_restart_hydro')
  115 FORMAT (' Equation of state tables have been loaded in setup_restart_hydro')
  117 FORMAT (' Composition mass fractions have been loaded for NSE material in setup_restart_hydro')
  119 FORMAT (' Convection variables have been initialized in setup_restart_hydro')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                 \\\\\ INITIALIZE VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize shifted radial array index boundaries
!-----------------------------------------------------------------------

ncycle                       = nrst
jr_min                       = imin + 1
jr_max                       = imax + 1
jr_maxp                      = jr_max + 1

!-----------------------------------------------------------------------
!  Set quatities at inner edge of configuration
!-----------------------------------------------------------------------

r(1)                         = zero

!-----------------------------------------------------------------------
!  Compute pseudoviscous multipliers
!-----------------------------------------------------------------------

DO j    = jr_min,jr_max
  IF ( q0_x(j) < zero ) q0_x(j) = q0_x(1)
END DO

DO j    = 1,ny
  IF ( q0_y(j) < zero ) q0_y(j) = q0_y(1)
END DO

DO j    = 1,nz
  IF ( q0_z(j) < zero ) q0_z(j) = q0_z(1)
END DO

WRITE (nlog,101)

!-----------------------------------------------------------------------
!  Transfer zone-centered x-velocity to shifted MGFLD arrays
!-----------------------------------------------------------------------

u(imin+1:imax+1)             = up(imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Set m-1 and m+1 values of velocities and set prescribed
!   values (if option selected).
!-----------------------------------------------------------------------

IF ( iuset == 1 ) THEN
  u(jr_min:jr_max)           = uset(jr_min:jr_max)
END IF

IF ( iubcjmx == 1 ) THEN
  u(jr_max)                  = ubcjmx
END IF

IF ( iubcjmn == 1 ) THEN
  u(1)                       = ubcjmn
END IF

WRITE (nlog,103)

!-----------------------------------------------------------------------
!  Transfer zone-centered x-velocity back to radial_array_module
!-----------------------------------------------------------------------

 up(imin:imax,ij_ray,ik_ray) = u(imin+1:imax+1)

!-----------------------------------------------------------------------
!
!                  \\\\\ TRANSFER VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer zone-centered independent variables to shifted mgfld arrays
!-----------------------------------------------------------------------

rho   (imin+1:imax+1)        = rhop   (imin:imax,ij_ray,ik_ray)
t     (imin+1:imax+1)        = tp     (imin:imax,ij_ray,ik_ray)
ye    (imin+1:imax+1)        = yep    (imin:imax,ij_ray,ik_ray)
be_nuc_rep(imin+1:imax+1)    = be_nuc_repp(imin:imax,ij_ray,ik_ray)
a_nuc_rep (imin+1:imax+1)    = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep (imin+1:imax+1)    = z_nuc_repp (imin:imax,ij_ray,ik_ray)

xn(imin+1:imax+1,:)          = xnp(imin:imax,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer zone-edgeed independent variables to shifted mgfld arrays
!-----------------------------------------------------------------------

r(imin:imax+1)               = rp(imin:imax+1)
r(imax-imin+3)               = r(imax-imin+2) + ( r(imax-imin+2) - r(imax-imin+1) )

WRITE (nlog,105)

!-----------------------------------------------------------------------
!  Modify problem before eos table setup.
!  The modifications are contained in subroutine pblmst1. If no
!   modifications, subroutine pblmst1 is a dummy subroutine.
!-----------------------------------------------------------------------

CALL pblmst1( ij_ray, ik_ray )
WRITE (nlog,107)

!-----------------------------------------------------------------------
!
!               \\\\\ SET EQUATION OF STATE TABLES ////
!
!  Load equation of state tables and set m-1 and m+1 values of
!   independent variables.
!-----------------------------------------------------------------------

rho(jr_maxp)                 = rho(jr_max)
t  (jr_maxp)                 = t(jr_max)
ye (jr_maxp)                 = ye(jr_max)

rhor(jr_min:jr_maxp)         = rho(jr_min:jr_maxp)
tr  (jr_min:jr_maxp)         = t(jr_min:jr_maxp)
yer (jr_min:jr_maxp)         = ye(jr_min:jr_maxp)

WRITE (nlog,109)
CALL esrgnz_x( jr_min, jr_maxp, rho, t, ye, ij_ray,ik_ray)
WRITE (nlog,111)
CALL eqstz_x ( jr_min, jr_maxp, rho, t, ye, ij_ray,ik_ray)
WRITE (nlog,113)
CALL gammaz_x( jr_min, jr_max , rho, t, ij_ray,ik_ray)

IF ( rho0(2) == zero ) THEN
  rho0(jr_min:jr_maxp)       = rho (jr_min:jr_maxp)
  ye0 (jr_min:jr_maxp)       = ye  (jr_min:jr_maxp)
  p0  (jr_min:jr_maxp)       = aesv(jr_min:jr_maxp,1,ij_ray,ik_ray)
END IF

WRITE (nlog,109)

!-----------------------------------------------------------------------
!  Put nse composition in composition arrays
!-----------------------------------------------------------------------

IF ( first_c ) THEN
  first_c                    = .false.
  i_n                        = nuc_number + 1
  i_p                        = nuc_number + 1
  i_4He                      = nuc_number + 1
  DO i = 1,nuc_number
    IF ( a_name(i) == '  n  ' ) THEN
      i_n                    = i
    END IF ! a_name(i) == '  n  '
    IF ( a_name(i) == '  p  ' ) THEN
      i_p                    = i
    END IF ! a_name(i) == '  p  '
    IF ( a_name(i) == '  4He' ) THEN
      i_4He                  = i
    END IF ! a_name(i) == '  4He'
  END DO ! i = 1,nuc_number
END IF ! first_c

DO j = jr_min,jr_max
  IF ( nse(j,ij_ray,ik_ray) == 0 ) CYCLE
  xn(j,:)                    = zero
  xn(j,i_n)                  = aesv(j,7,ij_ray,ik_ray)
  xn(j,i_p)                  = aesv(j,8,ij_ray,ik_ray)
  xn(j,nuc_number+1)         = aesv(j,9,ij_ray,ik_ray)
  xn(j,i_4He)                = DMAX1( one - xn(j,i_n) - xn(j,i_p) - xn(j,nuc_number+1), zero )
END DO ! j = jr_min,jr_max
WRITE (nlog,117)

!-----------------------------------------------------------------------
!  Transfer interpolated EOS quantities for storage in 
!-----------------------------------------------------------------------

aesv_c(imin:imax+1,:,ij_ray,ik_ray) = aesv(jr_min:jr_maxp,:,ij_ray,ik_ray)
xnp(imin:imax,:,ij_ray,ik_ray)      = xn(jr_min:jr_max,:)

!-----------------------------------------------------------------------
!  Pressure boundary condition
!-----------------------------------------------------------------------

IF ( jr_max /= jmr ) THEN
  IF ( ipbnd    == 1 ) aesv(jr_max+1,1,ij_ray,ik_ray) = aesv(jr_max,1,ij_ray,ik_ray)       &
&                            * aesv(jr_max,1,ij_ray,ik_ray)/aesv(jr_max-1,1,ij_ray,ik_ray)
  pbound                     = aesv(jr_max+1,1,ij_ray,ik_ray)
END IF

IF ( ipbnd    == 2 ) THEN
  aesv(jr_max+1,1,ij_ray,ik_ray) = zero
  pbound                     = zero
END IF

!-----------------------------------------------------------------------
!  Reset jmr
!-----------------------------------------------------------------------

jmr                          = jr_max

!-----------------------------------------------------------------------
!  Initialize convection variables
!-----------------------------------------------------------------------

aledoux(jr_min:jr_max)       = zero
psclht (jr_min:jr_max)       = zero
lmix   (jr_min:jr_max)       = zero
ulmixl (jr_min:jr_max)       = zero
ulcnvct(jr_min:jr_max)       = zero
pcnvct (jr_min:jr_max)       = zero
scnvct (jr_min:jr_max)       = zero

WRITE (nlog,119)

RETURN
END SUBROUTINE setup_restart_hydro

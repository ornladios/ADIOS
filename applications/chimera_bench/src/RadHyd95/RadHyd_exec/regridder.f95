SUBROUTINE regridder( imin, imax, nx, ij_ray_dim, ny, ik_ray_dim, nz, nez, &
& nnu, nnc, nse_c, x_ef, dx_cf, x_cf, rhobar0, rho_c, t_c, ye_c, u_c, u_e, &
& v_c, w_c, psi0_c, regrid, rho_regrid, time, t_bounce, l_regrid )
!-----------------------------------------------------------------------
!
!    File:         regridder
!    Module:       regridder
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/18/07
!
!    Purpose:
!      To regrid the model data into imax mass zones of distribution 
!       determined by the distribution keys. Regridding is carried out
!       to third order in space.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!  paraset        : computes constants which are used to interpolate a parabola on a quantity
!  parabola       : computes parabolic fitting coefficients for a variable
!
!    Input arguments:
!  imin           : minimum unshifted radial grid index
!  imax           : maximum unshifted radial grid index
!  nx             : radial array extent
!  ij_ray_dim     : the number of y-zones on a processor before swapping with y
!  ny             : radial array extent
!  ik_ray_dim     : the number of z-zones on a processor before swapping with z
!  nz             : z-array extent
!  nez            : neutrino energy array extent
!  nnu            : neutrino flavor array extent
!  nnc            : composition array extent
!  nse_c          : nuclear statistical equilibrium flag
!  x_ef           : x grid zone face locations at cycle end [cm]
!  dx_cf          : x_ef(i+1) - x_ef(i) [cm]
!  x_cf           : x grid zone midpoint locations at cycle end [cm]
!  rhobar0        : angularly averaged density [g cm^{-3}]
!  rho_c          : density [g cm^{-3}]
!  t_c            : temperature: zone average [K]
!  ye_c           : electron fraction: zone averaged
!  u_c            : zone centered average velocity x direction [cm s^{-1}]
!  u_e            : x grid edge velocity for moving grid option [cm s^{-1}]
!  v_c            : zone centered average velocity y direction [cm s^{-1}]
!  w_c            : zone centered average velocity z direction [cm s^{-1}]
!  psi0_c         : zero moment of the neurino occupation probability
!  regrid         : regrid option switch
!  time           : elapsed time [s]
!  t_bounce       : time of core bounce [s]
!
!    Output arguments:
!  x_ef           : x grid zone face locations at cycle end [cm]
!  dx_cf          : x grid zone face locations at cycle end [cm]
!  x_cf           : x grid zone face locations at cycle end [cm]
!  rhobar0        : angularly averaged density [g cm^{-3}]
!  rho_c          : density [g cm^{-3}]
!  t_c            : temperature: zone average [K]
!  ye_c           : electron fraction: zone averaged
!  u_c            : zone centered average velocity x direction [cm s^{-1}]
!  u_e            : x grid edge velocity for moving grid option [cm s^{-1}]
!  v_c            : zone centered average velocity y direction [cm s^{-1}]
!  w_c            : zone centered average velocity z direction [cm s^{-1}]
!  psi0_c         : zero moment of the neurino occupation probability
!  dtnph          : current time step [s]
!  dtnmh          : preceding time step [s]
!  l_regrid       : regrid flag
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  cycle_module, edit_module, mgfld_remap_module, nu_energy_grid_module,
!  parallel_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero, half, one, third, frpith, frpi
USE physcnst_module, ONLY: pi, msolar

USE cycle_module, ONLY: ncycle
USE edit_module, ONLY: nprint, nlog
USE mgfld_remap_module, ONLY : psi0_re
USE nu_energy_grid_module, ONLY : nnugp
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER (len=2)                    :: regrid              ! regrid option switch

INTEGER, INTENT(in)                  :: imin                ! minimum unshifted radial grid index
INTEGER, INTENT(in)                  :: imax                ! maximum unshifted radial grid index

INTEGER, INTENT(in)                  :: nx                  ! radial zone extent
INTEGER, INTENT(in)                  :: ij_ray_dim          ! the number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)                  :: ny                  ! y-array extent
INTEGER, INTENT(in)                  :: ik_ray_dim          ! the number of z-zones on a processor before swapping with y
INTEGER, INTENT(in)                  :: nz                  ! z-array extent
INTEGER, INTENT(in)                  :: nez                 ! neutrino energy array extent
INTEGER, INTENT(in)                  :: nnu                 ! neutrino flavor array extent
INTEGER, INTENT(in)                  :: nnc                 ! composition array extent

INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: nse_c      ! nuclear statistical equilibrium flag

REAL(KIND=double), INTENT(in)                             :: rho_regrid ! regrid up to density rho_regrid or 5 zones behind shock, whichever is larger [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx)              :: rhobar0    ! angularly averaged density [g cm^{-3}]
REAL(KIND=double), INTENT(in)                             :: time       ! elapsed time [s]
REAL(KIND=double), INTENT(in)                             :: t_bounce   ! time of core bounce [s]

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

LOGICAL, INTENT(out)                          :: l_regrid   ! regrid flag

REAL(KIND=double), INTENT(out), DIMENSION(nx) :: dx_cf      ! radial grid widths [cm]
REAL(KIND=double), INTENT(out), DIMENSION(nx) :: x_cf       ! radial grid zone midpoints

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx+1)                             :: x_ef     ! radial grid edges [cm]
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rho_c    ! zone densities of the rezoned model [g cm^{-3}]
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: t_c      ! temperature: zone average [K]
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: ye_c     ! electron fraction: zone averaged
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: u_c      ! zone centered average velocity x direction [cm s^{-1}]
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: u_e      ! x grid edge velocity for moving grid option [cm s^{-1}]
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: v_c      ! zone centered average velocity y direction [cm s^{-1}]
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: w_c      ! zone centered average velocity z direction [cm s^{-1}]
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0_c   ! zone centered average velocity z direction [cm s^{-1}]

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=14)                        :: var_name

LOGICAL                                   :: first = .true.

INTEGER                                   :: istat            ! allocation status

INTEGER                                   :: it               ! iteration index
INTEGER                                   :: it_count         ! iteration counter
INTEGER, PARAMETER                        :: it_max = 50      ! maximum number of iteration attempts
INTEGER                                   :: i_nonconv        ! number of consecutive failure to converge in density

INTEGER, PARAMETER                        :: ngeom = 0        ! geometry glag

INTEGER                                   :: i                ! radial zone index
INTEGER                                   :: ii               ! radial zone index
INTEGER                                   :: ij_ray           ! j-index of a radial ray
INTEGER                                   :: ik_ray           ! k-index of a radial ray
INTEGER                                   :: k                ! neutrino energy index
INTEGER                                   :: n                ! padded index array
INTEGER                                   :: nn               ! neutrino flavor index
INTEGER                                   :: nmin             ! minimum padded radial array index
INTEGER                                   :: nmax             ! maximum padded radial array index
INTEGER                                   :: ntot             ! nmax + 6

INTEGER                                   :: i_out            ! outermost radial zone index of region to be rezoned
INTEGER                                   :: n_out            ! shifted outermost radial zone index of region to be rezoned

INTEGER                                   :: i_10_hi          ! exponent of the upper density decade in the region to be re-gridded
INTEGER                                   :: i_10_lo          ! exponent of the lower density decade in the region to be re-gridded
INTEGER                                   :: n_scale          ! number of density scale heights in the region to be re-gridded
INTEGER, DIMENSION(100)                   :: i_bounds         ! indices of the density interval boundaries

INTEGER                                   :: ne               ! X_old(ne) > X_new(i)
INTEGER                                   :: nem1             ! ne - 1
INTEGER                                   :: nep              ! do index
INTEGER                                   :: ne_min           ! minimum index in the new grid to search for the next zone edge
INTEGER                                   :: nminn1           ! nmin+n-1

INTEGER                                   :: jjshockmx        ! maximum radial zone index of shock for a given ray
INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim) :: jjshockmn        ! minimum radial zone index of shock for a given ray
INTEGER, DIMENSION(ny,nz)                 :: jjshockmn_a      ! minimum radial zone index of shock for a given ray
INTEGER                                   :: ishockmin        ! minimum radial zone index
INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim) :: inse_mn          ! outermost radial zone index for which nes = 1
INTEGER, DIMENSION(ny,nz)                 :: inse_mn_a        ! outermost radial zone index for which nes = 1
INTEGER                                   :: inse_min         ! minimum radial zone index

REAL(KIND=double), PARAMETER                     :: tol = 1.d-6     ! iteration tolerence
REAL(KIND=double), PARAMETER                     :: pqmin = 1.d0    ! minimum ratio of pseudoviscosity to pressure to locate shock

REAL(KIND=double)                                :: thpifrth        ! 3/4pi

REAL(KIND=double)                                :: rho_ratio       ! final density ratio between adjacent sones in the regridded region
REAL(KIND=double)                                :: rho_min         ! regrid up to density rho_min [g cm^{-3}]
REAL(KIND=double)                                :: rhobar0_shk     ! mean density of inner shock radius
REAL(KIND=double)                                :: rho_ratio_b     ! ratio of rho(nmax) to rho(nmax-1) for computing outer ghosts
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rhobar_rgrd     ! mean densities to be achieved in the regridded region [g cm^{-3}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rhobar          ! mean densities actually achieved in the regridded region [g cm^{-3}]

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: mbar            ! enclosed mass of the rezoned model [g]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dmbar           ! zone mass of the rezoned model [g]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: m               ! enclosed mass of the rezoned model for a particular radial ray [g]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dm              ! zone mass of the rezoned modell for a particular radial ray [g]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: r0              ! zone radii of the model before rezoning [cm]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: r               ! zone radii of the rezoned model [cm]
REAL(KIND=double)                                :: r_ratio         ! ratio of r(i+1) to r(i) in regions where rho convergence fails
REAL(KIND=double), PARAMETER                     :: dr_min = 1.d4   ! minimum allowable r(i+1) - r(i) [cm]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: vol             ! rezoned enclosed volume [cm^{3}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol            ! rezoned zone volume [cm^{3}]

REAL(KIND=double)                                :: dvolume_min     ! minimum volume for bisection iteration
REAL(KIND=double)                                :: dvolume_max     ! maximum volume for bisection iteration
REAL(KIND=double)                                :: dvolume_test    ! test volume for bisection iteration
REAL(KIND=double)                                :: dvolume_2_nem1  ! volume between new zonez 2 and old zone nem1
REAL(KIND=double)                                :: dvolume_i_nem1  ! volume between new zonez i and pld zone nem1
REAL(KIND=double)                                :: dvolume0        ! volume of one or multiple old shells
REAL(KIND=double)                                :: d_vol_i_ne_min  ! volume between old zone ne_min and new zone i

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: vol_0           ! original enclosed volume [cn^{3}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dvol_0          ! original zone volume [cn^{3}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: xaV0            ! padded array of old enclosed volumes
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dxV0            ! padded array of old zone volumes
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: paraV           ! parabolic coefficients for volumes
REAL(KIND=double)                                :: V_ratio_b       ! ratio of dxV0(nmax) to dxV0(nmax-1) for computing outer ghosts

REAL(KIND=double)                                :: fourthd         ! 4/3
REAL(KIND=double)                                :: fractn          ! half the zone Lagrangian and Eulerian frid overlap
REAL(KIND=double)                                :: fractn2         ! 1 - 4/3*fractn
REAL(KIND=double)                                :: fracti          ! fraction of overlap to begin average
REAL(KIND=double)                                :: fractim1        ! fraction of overlap to begin average

REAL(KIND=double)                                :: dmass_2_nem1    ! mass between new zone 2 and old zone nem1
REAL(KIND=double)                                :: dmass_i_nem1    ! mass between new zone i and old zone nem1
REAL(KIND=double)                                :: dm_i_ne_min     ! mass between new zone i and old zone ne_min
REAL(KIND=double)                                :: dmass0          ! mass of one or multiple old shells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: m0              ! unshifted array of old enclosed masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dm0             ! unshifted array of old shell masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: xam0            ! padded array of old enclosed masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dxm0            ! padded array of old zone masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: param           ! parabolic coefficients for masses

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rho0            ! densities of the model before rezoning [g cm^{-3}] for a given ray
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rho             ! densities of the model after rezoning [g cm^{-3}] for a given ray
REAL(KIND=double)                                :: rho_test        ! dmass0/dvolume0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rho_pad         ! padded density array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rhol            ! density at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rho6            ! density parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: drho            ! density slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t0              ! temperature of the model before rezoning [K] for a given ray
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t               ! temperature of the model after rezoning [K] for a given ray
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t_pad           ! padded temperature array [K]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: tl              ! temperature at left edge of zone [K]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: t6              ! temperature parabola coeffecient [K]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dt              ! temperature slope [K]
REAL(KIND=double)                                :: t_lt            ! left part of average temperature [K]
REAL(KIND=double)                                :: t_rt            ! right part of average temperature [K]
REAL(KIND=double)                                :: t_incl          ! temperature between xaV0(ne_m) and xaV0(ne-1) [K]

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye0             ! electron fraction of the model before rezoning for a given ray
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye              ! electron fraction of the model after rezoning for a given ray
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye_pad          ! padded electron fraction array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: yel             ! electron fraction at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ye6             ! electron fraction parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dye             ! electron fraction slope
REAL(KIND=double)                                :: ye_lt           ! left part of average electron fraction
REAL(KIND=double)                                :: ye_rt           ! right part of average electron fraction
REAL(KIND=double)                                :: ye_incl         ! electron number between xaV0(ne_m) and xaV0(ne-1)

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: u0              ! x-velocity of the model before rezoning [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: u               ! x-velocity of the model after rezoning [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: u_pad           ! padded x-velocity array [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: ul              ! x-velocity at left edge of zone [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: u6              ! x-velocity parabola coeffecient [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: du              ! x-velocity slope [cm s^{-1}]
REAL(KIND=double)                                :: u_lt            ! left part of average x-velocity [cm s^{-1}]
REAL(KIND=double)                                :: u_rt            ! right part of average x-velocity [cm s^{-1}]
REAL(KIND=double)                                :: u_incl          ! x-velocity between xaV0(ne_m) and xaV0(ne-1) [cm s^{-1}]

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: v0              ! y-velocity of the model before rezoning [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: v               ! y-velocity of the model after rezoning [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: v_pad           ! padded y-velocity array [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: vl              ! y-velocity at left edge of zone [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: v6              ! y-velocity parabola coeffecient [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dv              ! y-velocity slope [cm s^{-1}]
REAL(KIND=double)                                :: v_lt            ! left part of average y-velocity [cm s^{-1}]
REAL(KIND=double)                                :: v_rt            ! right part of average y-velocity [cm s^{-1}]
REAL(KIND=double)                                :: v_incl          ! y-velocity between xaV0(ne_m) and xaV0(ne-1) [cm s^{-1}]

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: w0              ! z-velocity of the model before rezoning [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: w               ! z-velocity of the model after rezoning [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: w_pad           ! padded z-velocity array [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: wl              ! z-velocity at left edge of zone [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: w6              ! z-velocity parabola coeffecient [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dw              ! z-velocity slope [cm s^{-1}]
REAL(KIND=double)                                :: w_lt            ! left part of average z-velocity [cm s^{-1}]
REAL(KIND=double)                                :: w_rt            ! right part of average z-velocity [cm s^{-1}]
REAL(KIND=double)                                :: w_incl          ! z-velocity between xaV0(ne_m) and xaV0(ne-1) [cm s^{-1}]

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi0_0          ! zero neutrino moment before rezoning
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi0            ! zero neutrino moment after rezoning
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psil            ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi6            ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dpsi            ! zero neutrino moment slope
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: psi             ! 1-D zero neutrino moment array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psil_a          ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: psi6_a          ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:) :: dpsi_a          ! izero neutrino moment slope
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: psi_lt          ! left part of average z-velocity [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: psi_rt          ! right part of average z-velocity [cm s^{-1}]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: psi_incl        ! z-velocity between xaV0(ne_m) and xaV0(ne-1) [cm s^{-1}]

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dum             ! dummy array

REAL(KIND=double)                                :: rho_lft         ! left part of average density
REAL(KIND=double)                                :: rho_rgh         ! right part of average density

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Cannot find i_out in subroutine regridder; rhobar0(1)=',es11.3,' rhobar0(imax)=',es11.3,' rho_min=',es11.3)
  103 FORMAT (' nse(',i4,',',i4,',',i4,')=',i4,' cannot regrid non-NSE regions')
  111 FORMAT (' Iteration will not converge for rho_test in subroutine regrid for determining zone ',i4,/ &
&      ' dvolume_test=',es11.3,' dvol_0(nem1)=',es11.3,' rho_test=',es11.3,' rhobar_rgrd(i)=',es11.3)
  113 FORMAT (' Iteration will not converge for rho_test in subroutine regrid for determining zone ',i4,/ &
&      ' for the case in which i and i-1 are in the same original zone',/ &
&      ' dvolume_test=',es11.3,' dvol_0(nem1)=',es11.3,' rho_test=',es11.3,' rhobar_rgrd(i)=',es11.3)
  121 FORMAT (' Regridding is being attempted in regridder from the core center to rho_min=',es11.3, &
& ' ishockmin=',i4,' inse_min=',i4,' i_out=',i4,' rhobar0_shk=',es11.3,' rho_regrid=',es11.3) 
  123 FORMAT (' Regrid performed at time', es15.8,' s, time from bounce',es15.8,' cycle',i8)
 1001 FORMAT (' Allocation problem for array ',a10,' in regrid')
 2001 FORMAT (' Deallocation problem for array ',a10,' in regrid')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if regrid = 'no'
!-----------------------------------------------------------------------

l_regrid                 = .false.
IF ( regrid == 'no' ) RETURN
l_regrid                 = .true.

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (rhobar_rgrd(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhobar_rgrd   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rhobar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhobar        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (mbar(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'mbar          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmbar(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmbar         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (m(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'm             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dm(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r0            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r0(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (vol(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vol           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol          '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (vol_0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vol_0         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol_0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol_0        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xaV0(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xaV0          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dxV0(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dxV0          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (paraV(10,nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paraV         '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (m0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'm0            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dm0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xam0(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xam0          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dxm0(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dxm0          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (param(10,nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'param         '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rho0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho0          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rho(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rho_pad(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_pad       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rhol(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhol          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rho6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho6          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (drho(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'drho          '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (t0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't0            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t_pad(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_pad         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (tl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tl            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (t6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't6            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dt(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dt            '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ye0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye0           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye_pad(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_pad        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (yel(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yel           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye6           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dye(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye           '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (u0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u0            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u_pad(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_pad         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ul(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ul            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u6            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (du(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'du            '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (v0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v0            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v_pad(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_pad         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (vl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vl            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v6            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dv(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dv            '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (w0(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w0            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w_pad(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_pad         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (wl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wl            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w6            '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dw(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dw            '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi0_0(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_0        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi0(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psil(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psil_a(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil_a        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi6_a(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6_a        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi_a(nx+12,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi_a        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi_lt(nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi_lt        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi_rt(nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi_rt        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi_incl(nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi_incl      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dum(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dum           '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ DETERMINE DENSITY OUT TO REGRID /////
!
!-----------------------------------------------------------------------

DO ij_ray = 1, ij_ray_dim
  DO ik_ray = 1, ik_ray_dim
    CALL findshock_min( imin+1, imax+1, ij_ray, ik_ray, pqmin, &
&    jjshockmn(ij_ray,ik_ray), jjshockmx )
    inse_mn(ij_ray,ik_ray) = MINLOC( nse_c(:,ij_ray,ik_ray), DIM = 1 ) - 1
  END DO ! ik_ray = 1, ik_ray_dim
END DO ! ij_ray = 1, ij_ray_dim

jjshockmn_a(:,:)           = jjshockmn(:,:)
inse_mn_a(:,:)             = inse_mn(:,:)
ishockmin                  = MINVAL( jjshockmn(:,:) ) - 1
inse_min                   = MINVAL( inse_mn_a(:,:) )
i_out                      = MIN( ishockmin, inse_min )
rhobar0_shk                = rhobar0(i_out)
rho_min                    = DMAX1( rhobar0_shk, rho_regrid )

WRITE (nlog,121) rho_min, ishockmin, inse_min, i_out, rhobar0_shk, rho_regrid

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

fourthd                    = 4.d0/3.d0
thpifrth                   = 1.d0/frpith

rhobar_rgrd                = zero
rhobar                     = zero

mbar                       = zero
dmbar                      = zero
m                          = zero
dm                         = zero
r                          = zero
vol                        = zero
dvol                       = zero

r0                         = zero
vol_0                      = zero
dvol_0                     = zero
xaV0                       = zero
dxV0                       = zero
paraV                      = zero

m0                         = zero
dm0                        = zero
xam0                       = zero
dxm0                       = zero
param                      = zero

rho0                       = zero
rho                        = zero
rho_pad                    = zero
rhol                       = zero
rho6                       = zero
drho                       = zero

t0                         = zero
t                          = zero
t_pad                      = zero
tl                         = zero
t6                         = zero
dt                         = zero

ye0                        = zero
ye                         = zero
ye_pad                     = zero
yel                        = zero
ye6                        = zero
dye                        = zero

u0                         = zero
u                          = zero
u_pad                      = zero
ul                         = zero
u6                         = zero
du                         = zero

v0                         = zero
v                          = zero
v_pad                      = zero
vl                         = zero
v6                         = zero
dv                         = zero

w0                         = zero
w                          = zero
w_pad                      = zero
wl                         = zero
w6                         = zero
dw                         = zero

psi0_0                     = zero
psi0                       = zero
psil                       = zero
psi6                       = zero
dpsi                       = zero
psil_a                     = zero
psi6_a                     = zero
dpsi_a                     = zero
psi_lt                     = zero
psi_rt                     = zero
psi_incl                   = zero

dum                        = zero

r0(:)                      = x_ef(:)

!-----------------------------------------------------------------------
!  Calculate original zone volumes
!-----------------------------------------------------------------------

vol_0(1)                   = zero
vol_0(2:imax+1)            = frpith * x_ef(2:imax+1)**3
dvol_0(1:imax)             = vol_0(2:imax+1) - vol_0(1:imax)

!-----------------------------------------------------------------------
!  Calculate original mass shells and enclosed masses
!-----------------------------------------------------------------------

dm0(1:imax)                = rhobar0(1:imax) * dvol_0(1:imax)

m0(1)                      = zero
DO i = 1,imax
  m0(i+1)                  = m0(i) + dm0(i)
END DO ! i = 1,imax

!-----------------------------------------------------------------------
!  Determine i_out such that rhobar0(i_out) < rho_min < rhobar0(i_out+1)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                  |                     |                   |
!                  |                     |                   |
!                  |                     |<-----rhobar0----->|
!                  |                     |                   |
!                  |                     |                   |
!
!               i_out-1                i_out             i_out+1
!
!-----------------------------------------------------------------------

DO i = 1, imax
  IF ( rhobar0(i) >= rho_min  .and.  rhobar0(i+1) <= rho_min ) THEN
    i_out                  = i
    EXIT
  END IF ! rhobar0(i) >= rho_min  .and.  rhobar0(i+1) <= rho_min
END DO ! i = 1, imax

IF ( i_out == imax ) THEN
  WRITE (nlog,101) rhobar0(1), rhobar0(imax), rho_min
  STOP
END IF ! i_out == imax

DO ij_ray = 1, ij_ray_dim
  DO ik_ray = 1, ik_ray_dim
    IF ( nse_c(i_out,ij_ray,ik_ray) == 0 ) THEN
      WRITE (nlog,103) i_out, ij_ray, ik_ray, nse_c(i_out,ij_ray,ik_ray)
      GO TO 5000
    END IF ! nse_c(i_out,ij_ray,ik_ray) == 1
  END DO ! ik_ray = 1, ik_ray_dim
END DO ! ij_ray = 1, ij_ray_dim

nmin                       = 7
nmax                       = imax + 6
ntot                       = nmax + 6
n_out                      = i_out + 6

!-----------------------------------------------------------------------
!  rho_ratio is the ratio between adjacent zones to be achieved in the
!   regridded region
!
!              i_out - 1
!     rho_ratio = rhobar0(i_out)/rhobar0(1)
!-----------------------------------------------------------------------

rho_ratio                  = ( rhobar0(i_out)/rhobar0(1) )**( 1.d0/DBLE(i_out-1) )

!-----------------------------------------------------------------------
!  Now calculate the mean densities in the regridded region
!-----------------------------------------------------------------------

rhobar_rgrd(1)             = rhobar0(1)
DO i = 2, i_out
  rhobar_rgrd(i)           =  rhobar_rgrd(i-1) * rho_ratio
END DO

!-----------------------------------------------------------------------
!
!         \\\\\ PPM INTERPOLATE IN VOLUME THE DEMSITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Put original volume variables into arrays, padding with 6 ghost zones
!        
!  dxV(n) must be indexed as the shell above enxlosed volume xaV(n)
!-----------------------------------------------------------------------

xaV0(nmin:nmax+1)          = vol_0(imin:imax+1)
dxV0(nmin:nmax)            = dvol_0(imin:imax)

!-----------------------------------------------------------------------
!  Fill ghost zones
!-----------------------------------------------------------------------

V_ratio_b                  = dxV0(nmax)/dxV0(nmax-1)
DO n = 1, 6
  dxV0(nmin-n)             = dxV0(nmin+n-1)
  xaV0(nmin-n)             = xaV0(nmin-n+1) - dxV0(nmin-n)
  dxV0(nmax+n)             = dxV0(nmax) * V_ratio_b**(n)
  xaV0(nmax+n)             = xaV0(nmax+n-1) + dxV0(nmax+n)
END DO ! n = 1, 6

!-----------------------------------------------------------------------
!  Compute parabolic coefficients in volume
!-----------------------------------------------------------------------

CALL paraset( nx+12, paraV, dxV0, xaV0, nmin-4, nmax+4, ngeom )

!-----------------------------------------------------------------------
!  Put original mean densities into arrays, padding with 6 ghost zones
!
!  rho_pad(n), etc. must be indexed as the shells above enxlosed mass
!   xam(n)
!-----------------------------------------------------------------------

rho_pad(nmin:nmax)         = rhobar0(imin:imax)

!-----------------------------------------------------------------------
!  Fill ghost zones
!-----------------------------------------------------------------------

rho_ratio_b                = rho_pad(nmax)/rho_pad(nmax-1)
DO n = 1, 6
  rho_pad(nmin-n)          = rho_pad(nmin+n-1)
  rho_pad(nmax+n)          = rho_pad(nmax) * rho_ratio_b**(n)
END DO

!-----------------------------------------------------------------------
!  Generate parabolic interpolation coefficients
!-----------------------------------------------------------------------

CALL parabola( nmin-3, nmax+3, ntot, paraV, rho_pad, drho, rho6, rhol, &
& dum, 0, 0, ngeom )

!-----------------------------------------------------------------------
!
!         \\\\\ DETERMINE THE LOCATIONS OF THE NEW ZONES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                   ||||| Radial zone 1 |||||
!-----------------------------------------------------------------------

rhobar(1)                  = rhobar0(1)
mbar(1)                    = zero
dmbar(1)                   = dm0(1)
mbar(2)                    = dmbar(1)
r(1)                       = zero
r(2)                       = r0(2)
dvol(1)                    = dvol_0(1)
vol(1)                     = zero
vol(2)                     = frpith * r(2)**3

!-----------------------------------------------------------------------
!                   ||||| Radial zone 2 |||||
!-----------------------------------------------------------------------

ne_min                     = 2
i                          = 2

!-----------------------------------------------------------------------
!  Find position of new zone edge relative to the nearby zone edges
!   in the original grid.
!
!  ne is the index of the original grid for which the mean density
!   between i and ne falls below rhobar_rgrd(i)
!-----------------------------------------------------------------------

dmass_2_nem1               = zero
dvolume_2_nem1             = zero
DO nep = ne_min, i_out
  dmass_2_nem1             = dmass_2_nem1 + dm0(nep)
  dvolume_2_nem1           = dvolume_2_nem1 + dvol_0(nep)
  rho_test                 = dmass_2_nem1/dvolume_2_nem1
  IF ( rho_test < rhobar_rgrd(i) ) THEN
    ne                     = nep + 1
    nem1                   = nep
    dmass_2_nem1           = dmass_2_nem1 - dm0(nep)
    dvolume_2_nem1         = dvolume_2_nem1 - dvol_0(nep)
    EXIT
  END IF
END DO ! nep = ne_min, i_out

!-----------------------------------------------------------------------
!  Iterate to find the fraction of zone mem1 --> ne for which
!   rho_test = rhobar_rgrd(i)
!
!        :             :     |    :
!        :             :     |    :
!        :     ...     :     |    :
!        :             :     |    :
!        :             :     |    :
!
!    ne_min=2         nem1  i+1   ne
!        i
!                      <---- n --->
!
!-----------------------------------------------------------------------

n                          = nem1 + 6
dvolume_min                = zero
dvolume_max                = dvol_0(nem1)
DO it = 1, it_max

  it_count                 = it
  dvolume_test             = half * ( dvolume_min + dvolume_max )
  fractn                   = half * dvolume_test/dvol_0(nem1)
  fractn2                  = 1.d0 - fourthd * fractn
  rho_rgh                  = rhol(n) + fractn * ( drho(n) + fractn2 * rho6(n) )
  rho_test                 = ( dmass_2_nem1 + rho_rgh * dvolume_test )/( dvolume_2_nem1 + dvolume_test )
    
  IF ( DABS( rho_test - rhobar_rgrd(i) )/rhobar_rgrd(i) < tol ) EXIT

  IF ( rho_test >= rhobar_rgrd(i) ) THEN
    dvolume_min            = dvolume_test
  ELSE ! rho_test < rhobar_rgrd(i)
    dvolume_max            = dvolume_test
  END IF ! rho_test >= rhobar_rgrd(i)

END DO ! it = 1, it_max

IF ( it_count == it_max ) THEN
  WRITE (nlog,111) i, dvolume_test, dvol_0(nem1), rho_test, rhobar_rgrd(i)
  IF ( myid == 0 ) WRITE (nprint,111) i, dvolume_test, dvol_0(nem1), rho_test, rhobar_rgrd(i)
  STOP
END IF ! it == it_max
  
!-----------------------------------------------------------------------
!  Success! Find new shell mass, enclosed mass, etc. of the new mass
!   shell
!-----------------------------------------------------------------------
  
rhobar(i)                  = rho_test
dmbar(i)                   = dmass_2_nem1 + rho_rgh * dvolume_test
mbar(i+1)                  = mbar(i) + dmbar(i)
r(i+1)                     = ( r(i)**3 + thpifrth * ( dvolume_2_nem1 + dvolume_test ) )**third
vol(i+1)                   = frpith * r(i+1)**3
dvol(i)                    = vol(i+1) - vol(i)

!-----------------------------------------------------------------------
!              ||||| Radial zones 3 < i < i_out - 1 |||||
!-----------------------------------------------------------------------

Outer_i_loop1: DO i = 3, i_out - 1

!-----------------------------------------------------------------------
!  ne_min -1 and ne_min are the original grid edges that are on the
!   immediate left and right, respectively, of new grid edge i
!
!  After each loop over i, ne lies just to the right of i
!-----------------------------------------------------------------------
  DO nep = 1, i_out
    IF ( r0(nep) >= r(i) ) THEN
      ne                   = nep
      nem1                 = nep - 1
      EXIT
    END IF !  r0(nep) >= r(ii)
  END DO ! nep = 1, i_out

  ne_min                   = ne
  dmass0                   = zero
  dvolume0                 = zero

!-----------------------------------------------------------------------
!  Does i+1 lie between ne_min-1 and ne_min?
!
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    : 
!
!   ne_min-1    i         i+1  ne_min
!
!        <---------- n --------->
!
!               dvolume_test
!               <---------->
!
!                 d_vol_i_ne_min
!               <--------------->
!
!  IF so, iterate between i and ne_min to locate i+1
!-----------------------------------------------------------------------

  n                        = ne_min - 1 + 6
  d_vol_i_ne_min           = vol_0(ne_min) - vol(i)
  fractn                   = half * d_vol_i_ne_min/dvol_0(ne_min-1)
  fractn2                  = 1.d0 - fourthd * fractn
  rho_lft                  = rhol(n) + drho(n) - fractn * ( drho(n) - fractn2 * rho6(n) )

!-----------------------------------------------------------------------
!  If rho_lft <= rhobar_rgrd(i), i+1 lies between ne_min-1 and ne_min
!   to the right of i
!-----------------------------------------------------------------------

  IF ( rho_lft <= rhobar_rgrd(i) ) THEN

    dvolume_min            = zero
    dvolume_max            = d_vol_i_ne_min
    DO it = 1, it_max

      it_count             = it
      dvolume_test         = half * ( dvolume_min + dvolume_max )
      fractim1             = ( vol(i)                - vol_0(ne_min-1) )/dvol_0(ne_min-1)
      fracti               = ( vol(i) + dvolume_test - vol_0(ne_min-1) )/dvol_0(ne_min-1)
      rho_test             = rhol(n) + half * ( fractim1 + fracti ) * ( drho(n) + rho6(n) ) &
&                          - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * rho6(n)
  
      IF ( DABS( rho_test - rhobar_rgrd(i) )/rhobar_rgrd(i) < tol ) EXIT

      IF ( rho_test >= rhobar_rgrd(i) ) THEN
        dvolume_min        = dvolume_test
      ELSE ! rho_test < rhobar_rgrd(i)
        dvolume_max        = dvolume_test
      END IF ! rho_test >= rhobar_rgrd(i)

    END DO ! it = 1, it_max

    IF ( it_count == it_max ) THEN
      i_nonconv            = i_nonconv + 1

!-----------------------------------------------------------------------
!  If rho_lft <= DSQRT( rhobar_rgrd(i) * rhobar_rgrd(i-1) ), temporarly
!   locate i+1 at ne_min
!-----------------------------------------------------------------------

      IF ( rho_lft <= DSQRT( rhobar_rgrd(i) * rhobar_rgrd(i-1) ) ) THEN
        rho_test           = rho_lft
        dvolume_test       = vol_0(ne_min) - vol(i)
      ELSE ! rho_lft > DSQRT( rhobar_rgrd(i) * rhobar_rgrd(i-1) )
        WRITE (nlog,113) i, dvolume_test, d_vol_i_ne_min, rho_test, rhobar_rgrd(i)
        IF ( myid == 0 ) WRITE (nprint,113) i, dvolume_test, d_vol_i_ne_min, rho_test, rhobar_rgrd(i)
        STOP
      END IF ! rho_lft <= DSQRT( rhobar_rgrd(i) * rhobar_rgrd(i-1) )
    END IF ! it_count == it_max
  
!-----------------------------------------------------------------------
!  Success! Find new shell mass, enclosed mass, etc. of the new mass
!   shell
!-----------------------------------------------------------------------
  
    rhobar(i)              = rho_test
    dmbar(i)               = rho_test * dvolume_test
    mbar(i+1)              = mbar(i) + dmbar(i)
    r(i+1)                 = ( r(i)**3 + thpifrth * dvolume_test )**third
    vol(i+1)               = frpith * r(i+1)**3
    dvol(i)                = vol(i+1) - vol(i)

  ELSE ! rho_lft > rhobar_rgrd(i)

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i < ne
!  Iterate to find the fraction of zone mem1 --> ne for which
!   rho_test = rhobar_rgrd(i)
!
!        :        |     :             :     |    :
!        :        |     :             :     |    :
!        :        |     :     ...     :     |    :
!        :        |     :             :     |    :
!        :        |     :             :     |    :
!
!     ne_min-1    i   ne_min          ne-1 i+1   ne
!
!                 <---dvolume_i_nem1-->
!                 <----dmass_i_nem1--->
!
!-----------------------------------------------------------------------
!  Find position of new zone edge relative to the nearby zone edges
!   in the original grid.
!
!  ne is the index of the original grid for which the mean density falls
!   below rho_new(i)
!
!  dmass_i_nem1 is initially the mass from i to ne_min
!  dvolume_i_nem1 is initially the volume from i to ne_min, calculated above
!-----------------------------------------------------------------------

    dvolume_i_nem1         = d_vol_i_ne_min
    dmass_i_nem1           = rho_lft * d_vol_i_ne_min

    DO nep = ne_min, i_out

      dmass_i_nem1         = dmass_i_nem1 + dm0(nep)
      dvolume_i_nem1       = dvolume_i_nem1 + dvol_0(nep)
      rho_test             = dmass_i_nem1/dvolume_i_nem1

      IF ( rho_test < rhobar_rgrd(i) ) THEN
        ne                 = nep + 1
        nem1               = nep
        dmass_i_nem1       = dmass_i_nem1 - dm0(nep)
        dvolume_i_nem1     = dvolume_i_nem1 - dvol_0(nep)
        EXIT
      END IF

    END DO ! nep = ne_min, i_out

    n                      = nem1 + 6
    dvolume_min            = zero
    dvolume_max            = dvol_0(nem1)
    DO it = 1, it_max

      it_count             = it
      dvolume_test         = half * ( dvolume_min + dvolume_max )
      fractn               = half * dvolume_test/dvol_0(nem1)
      fractn2              = 1.d0 - fourthd * fractn
      rho_rgh              = rhol(n) + fractn * ( drho(n) + fractn2 * rho6(n) )
      rho_test             = ( dmass_i_nem1 + rho_rgh * dvolume_test )/( dvolume_i_nem1 + dvolume_test )
  
      IF ( DABS( rho_test - rhobar_rgrd(i) )/rhobar_rgrd(i) < tol ) EXIT

      IF ( rho_test >= rhobar_rgrd(i) ) THEN
        dvolume_min        = dvolume_test
      ELSE ! rho_test < rhobar_rgrd(i)
        dvolume_max        = dvolume_test
      END IF ! rho_test >= rhobar_rgrd(i)

    END DO ! it = 1, it_max

    IF ( it_count == it_max ) THEN
      i_nonconv            = i_nonconv + 1
      WRITE (nlog,113) i, dvolume_test, dvolume_i_nem1, rho_test, rhobar_rgrd(i)
      IF ( myid == 0 ) WRITE (nprint,113) i, dvolume_test, dvolume_i_nem1, rho_test, rhobar_rgrd(i)
    END IF ! it == it_max
  
!-----------------------------------------------------------------------
!  Success! Find new shell mass, enclosed mass, etc. of the new mass
!   shell
!-----------------------------------------------------------------------
  
    rhobar(i)              = rho_test
    dmbar(i)               = dmass_i_nem1 + rho_rgh * dvolume_test
    mbar(i+1)              = mbar(i) + dmbar(i)
    r(i+1)                 = ( r(i)**3 + thpifrth * ( dvolume_i_nem1 + dvolume_test ) )**third
    vol(i+1)               = frpith * r(i+1)**3
    dvol(i)                = vol(i+1) - vol(i)

  END IF ! rho_lft <= rhobar_rgrd(i)

!-----------------------------------------------------------------------
!       ||||| Radial zones i - i_nonconv + 1 < ii < i+1 |||||
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!              i                  i+1      i+i_nonconv  i+i_nonconv+1
!
!            conv               nonconv       nonconv          conv
!    :         |        :          |             |              |    :
!    :         |        :          |             |              |    :
!    :         |        :          |     ...     |              |    :
!    :         |        :          |             |              |    :
!    :         |        :          |             |              |    :
!
! ne_min-1 i-i_nonconv  ne    i-i_nonconv+1      i             I+1   
!              ii                 ii+1
!
!  If prior rho convergence failures, space the radii geometrically in
!   the region where rho convergence failed.
!-----------------------------------------------------------------------

  IF ( it_count /= it_max  .and.  i_nonconv /= 0 ) THEN
    r_ratio                = ( r(i+1)/r(i-i_nonconv) )**( 1.d0/DBLE(i_nonconv+1) )
    DO ii = 1, i_nonconv
      r(i-i_nonconv+ii)    = r(i-i_nonconv+ii-1) * r_ratio
      vol(i-i_nonconv+ii)  = frpith * r(i-i_nonconv+ii)**3
      dvol(i-i_nonconv+ii-1) = vol(i-i_nonconv+ii) - vol(i-i_nonconv+ii-1)
    END DO ! ii = 1, i_nonconv
    dvol(i)                = vol(i+1) - vol(i)

!-----------------------------------------------------------------------
!  Now recompute the rho's and dm's
!-----------------------------------------------------------------------

    DO ii = i - i_nonconv, i

!-----------------------------------------------------------------------
!              i                  i+1      i+i_nonconv  i+i_nonconv+1
!
!            conv               nonconv       nonconv          conv
!    :         |        :          |             |              |    :
!    :         |        :          |             |              |    :
!    :         |        :          |     ...     |              |    :
!    :         |        :          |             |              |    :
!    :         |        :          |             |              |    :
!
! ne_min-1 i-i_nonconv  ne    i-i_nonconv+1      i             I+1   
!              ii                 ii+1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find ne such that r0(ne-1) < r(ii) < r0(ne)
!-----------------------------------------------------------------------

      DO nep = 1, i_out
        IF ( r0(nep) >= r(ii) ) THEN
          ne               = nep
          nem1             = nep - 1
          EXIT
        END IF !  r0(nep) >= r(ii)
      END DO ! nep = 1, i_out
          
      ne_min               = ne
      dmass0               = zero
      dvolume0             = zero

      IF ( r(ii+1) <= r0(ne_min) ) THEN

!-----------------------------------------------------------------------
!  Does ii+1 lie between ne_min-1 and ne_min?
!
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    : 
!
!   ne_min-1    ii       ii+1  ne_min
!
!        <---------- n --------->
!
!  If so, average between ii and ii+1
!-----------------------------------------------------------------------

        n                  = ne_min - 1 + 6
        fractim1           = ( vol(ii  ) - vol_0(ne_min-1) )/dvol_0(ne_min-1)
        fracti             = ( vol(ii+1) - vol_0(ne_min-1) )/dvol_0(ne_min-1)
        rhobar(ii)         = rhol(n) + half * ( fractim1 + fracti ) * ( drho(n) + rho6(n) ) &
&                          - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * rho6(n)
        dmbar(ii)          = rhobar(ii) * dvol(ii)
        mbar(ii+1)         = mbar(ii) + dmbar(ii)

      ELSE ! r(ii+1) > r0(ne_min)

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i+1 < ne
!
!        :        |     :          :             :     |    :
!        :        |     :          :             :     |    :
!        :        |     :          :     ...     :     |    :
!        :        |     :          :             :     |    :
!        :        |     :          :             :     |    :
!
!     ne_min-1    ii  ne_min   ne_min+1         ne-1 ii+1   ne
!
!                 <--------dvolume_i_nem1------->
!                 <---------dmass_i_nem1-------->
!
!-----------------------------------------------------------------------
!  dmass_i_nem1 is initially the mass from i to ne_min
!  dvolume_i_nem1 is initially the volume from i to ne_min, calculated above
!-----------------------------------------------------------------------

        n                  = ne_min - 1 + 6
        d_vol_i_ne_min     = vol_0(ne_min) - vol(ii)
        fractn             = half * d_vol_i_ne_min/dvol_0(ne_min-1)
        fractn2            = 1.d0 - fourthd * fractn
        rho_lft            = rhol(n) + drho(n) - fractn * ( drho(n) - fractn2 * rho6(n) )
        dvolume_i_nem1     = d_vol_i_ne_min
        dmass_i_nem1       = rho_lft * d_vol_i_ne_min

        DO nep = ne_min, i_out
          IF ( r(ii+1) >= r0(nep+1) ) THEN
            dmass_i_nem1   = dmass_i_nem1 + dm0(nep)
            dvolume_i_nem1 = dvolume_i_nem1 + dvol_0(nep)
          ELSE ! r(ii+1) < r0(nep+1)
            ne             = nep + 1
            nem1           = nep
            n              = nem1 + 6
            fractn         = half * ( vol(ii+1) - vol_0(nem1) )/dvol_0(nem1)
            fractn2        = 1.d0 - fourthd * fractn
            rho_rgh        = rhol(n) + fractn * ( drho(n) + fractn2 * rho6(n) )
            dvolume_i_nem1 = dvolume_i_nem1 + vol(ii+1) - vol_0(nem1)
            dmass_i_nem1   = dmass_i_nem1 + rho_rgh * ( vol(ii+1) - vol_0(nem1) )
            EXIT
          END IF ! r(ii+1) >= r0(nep+1)
        END DO ! nep = ne_min, i_out

        dmbar(ii)          = dmass_i_nem1
        mbar(ii+1)         = mbar(ii) + dmbar(ii)
        rhobar(ii)         = dmbar(ii)/dvol(ii)

      END IF ! r(ii+1) <= r0(ne_min)

    END DO ! ii = i - i_nonconv, i
    i_nonconv              = 0
  END IF ! it_count /= it_max  .and.  i_nonconv /= 0

END DO Outer_i_loop1

!-----------------------------------------------------------------------
!                   ||||| Radial zones i_out |||||
!-----------------------------------------------------------------------

dmbar(i_out)               = m0(i_out+1) - mbar(i_out)
vol(i_out+1)               = vol_0(i_out+1)
dvol(i_out)                = vol(i_out+1) - vol(i_out)
rhobar(i_out)              = dmbar(i_out)/dvol(i_out)
mbar(i_out+1)              = m0(i_out+1)
r(i_out+1)                 = r0(i_out+1)

!-----------------------------------------------------------------------
!              ||||| Radial zones i_out+1 to imax |||||
!-----------------------------------------------------------------------

DO i = i_out+1, imax
  dmbar(i)                 = dm0(i)
  mbar(i+1)                = m0(i+1)
  vol(i+1)                 = vol_0(i+1)
  dvol(i)                  = dvol_0(i)
  rhobar(i)                = rhobar0(i)
  r(i+1)                   = x_ef(i+1)
END DO ! i = i_out+1, imax

!-----------------------------------------------------------------------
!
!               \\\\\ KEEP R(I+1) - R(I) > DR_MIN /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Prevent too narrow a zone width
!-----------------------------------------------------------------------

DO i = 2, i_out - 1
  IF ( r(i+1) - r(i) < dr_min ) r(i+1) = r(i) + dr_min
END DO

DO i = 2, i_out - 1
  vol(i+1)                 = frpith * r(i+1) * r(i+1) * r(i+1)
  dvol(i)                  = vol(i+1) - vol(i)
END DO
dvol(i_out)                = vol(i_out+1) - vol(i_out)

!-----------------------------------------------------------------------
!
!        \\\\\ NEW GRID HAS BEEN PRELIMINARILY DETERMINED /////
!
!-----------------------------------------------------------------------

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!-----------------------------------------------------------------------
!
!                     \\\\\ SMOOTH THE GRID /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Re-grid to constant dr/r for each decade in density
!-----------------------------------------------------------------------

i_10_hi                    = EXPONENT( rhobar0(1) )
i_10_lo                    = EXPONENT( rho_min    )
n_scale                    = i_10_hi - i_10_lo
DO i = i_10_lo, i_10_hi
  ii                       = i_10_hi - i + 1
  i_bounds(ii)             = MAXLOC( rhobar(:), DIM = 1, MASK = EXPONENT(rhobar) == i )
END DO
i_bounds(1)                = MAX( i_bounds(1), 2 )
i_bounds( i_10_hi-i_10_lo+2) = i_out

DO i = 1, n_scale + 1
WRITE (nlog,3011) i, i_bounds(i)
 3011 FORMAT (' i=',i4,' i_bounds(i)=',i4)
END DO

IF ( i_bounds( i_10_hi-i_10_lo+2) - i_bounds( i_10_hi-i_10_lo+1) >= 3 ) n_scale = n_scale + 1
DO i = 1, n_scale
  r_ratio                  = ( r(i_bounds(i+1))/r(i_bounds(i)) )**( 1.d0/DBLE( i_bounds(i+1) - i_bounds(i) ) )
  WRITE (nlog,3013) r(i_bounds(i)), r(i_bounds(i+1))
 3013 FORMAT (' r(i_bounds(i))=',es13.4,' r(i_bounds(i+1))=',es11.3)
  DO ii = i_bounds(i) + 1, i_bounds(i+1)
    r(ii)                  = r(ii-1) * r_ratio
  WRITE (nlog,3012) ii, r(ii)
 3012 FORMAT (' ii=',i4,' r(ii)=',es12.4)
  END DO
END DO

DO i = 2, i_out - 1
  vol(i+1)                 = frpith * r(i+1) * r(i+1) * r(i+1)
  dvol(i)                  = vol(i+1) - vol(i)
END DO
dvol(i_out)                = vol(i_out+1) - vol(i_out)

!-----------------------------------------------------------------------
!
!               \\\\\ NEW GRID HAS BEEN DETERMINED /////
!
!-----------------------------------------------------------------------

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!-----------------------------------------------------------------------
!
!        \\\\\ DETERMINE RHO AND M FOR EACH RADIAL RAY /////
!
!-----------------------------------------------------------------------


DO ij_ray = 1, ij_ray_dim
  DO ik_ray = 1, ik_ray_dim

!-----------------------------------------------------------------------
!  Store original mass arrays for ray ij_ray, ik_ray
!-----------------------------------------------------------------------

    rho0(:)                = rho_c(:,ij_ray,ik_ray)
    dm0(1:imax)            = rho0(1:imax) * dvol_0(1:imax)
    m0(1)                  = zero
    DO i = 1,imax
      m0(i+1)              = m0(i) + dm0(i)
    END DO ! 1,imax

!-----------------------------------------------------------------------
!  Put original mean densities into arrays, padding with 6 ghost zones
!
!  rho_pad(n), etc. must be indexed as the shells above enxlosed mass
!   xam(n)
!-----------------------------------------------------------------------

    rho_pad(nmin:nmax)     = rho_c(imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Fill ghost zones
!-----------------------------------------------------------------------

    rho_ratio_b            = rho_pad(nmax)/rho_pad(nmax-1)
    DO n = 1, 6
      rho_pad(nmin-n)      = rho_pad(nmin+n-1)
      rho_pad(nmax+n)      = rho_pad(nmax) * rho_ratio_b**(n)
    END DO ! n = 1, 6

!-----------------------------------------------------------------------
!  Generate parabolic interpolation coefficients
!-----------------------------------------------------------------------

    CALL parabola( nmin-3, nmax+3, ntot, paraV, rho_pad, drho, rho6, rhol, &
&    dum, 0, 0, ngeom )

!-----------------------------------------------------------------------
!                   ||||| Radial zone 1 |||||
!-----------------------------------------------------------------------

    rho(1)                 = rho0(1)
    m(1)                   = zero
    dm(1)                  = rho_c(1,ij_ray,ik_ray) * dvol_0(1)
    m(2)                   = dm(1)

!-----------------------------------------------------------------------
!                   ||||| Radial zone 2 |||||
!-----------------------------------------------------------------------

    ne_min                 = 2
    i                      = 2

!-----------------------------------------------------------------------
!  Find position of new zone edge relative to the nearby zone edges
!   in the original grid.
!
!  ne is the index of the original grid for which
!    r(ne-1) < r(i+1) < r(ne)
!-----------------------------------------------------------------------

    dmass_2_nem1           = zero
    dvolume_2_nem1         = zero
    DO nep = ne_min, i_out
      ne                   = nep + 1
      nem1                 = nep
      IF ( r0(nep+1) > r(i+1) ) EXIT
      dmass_2_nem1         = dmass_2_nem1 + rho0(nem1) * dvol_0(nem1)
      dvolume_2_nem1       = dvolume_2_nem1 + dvol_0(nem1)
    END DO ! nep = ne_min, i_out

!-----------------------------------------------------------------------
!  Determine the mass and volume of material between nem1 and i
!
!        :             :     |    :
!        :             :     |    :
!        :     ...     :     |    :
!        :             :     |    :
!        :             :     |    :
!
!    ne_min=2         nem1  i+1   ne
!        i
!                      <---- n --->
!
!-----------------------------------------------------------------------

    n                      = nem1 + 6
    fractn                 = half * ( vol(i+1) - vol_0(nem1) )/dvol_0(nem1)
    fractn2                = 1.d0 - fourthd * fractn
    rho_rgh                = rhol(n) + fractn * ( drho(n) + fractn2 * rho6(n) )
    dm(i)                  = dmass_2_nem1 + rho_rgh * ( vol(i+1) - vol_0(nem1) )
    rho(i)                 = dm(i)/dvol(i)
    m(i+1)                 = m(i) + dm(i)

!-----------------------------------------------------------------------
!              ||||| Radial zones 3 < i < i_out - 1 |||||
!-----------------------------------------------------------------------

    Outer_i_loop2: DO i = 3, i_out - 1

      ne_min               = ne
      dmass0               = zero
      dvolume0             = zero

!-----------------------------------------------------------------------
!  Does i+1 lie between ne_min-1 and ne_min?
!
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    : 
!
!   ne_min-1    i         i+1  ne_min
!
!        <---------- n --------->
!
!  If so, average between i and i+1
!-----------------------------------------------------------------------

      IF ( r(i+1) <= r0(ne_min) ) THEN
        n                  = ne_min - 1 + 6
        fractim1           = ( vol(i  ) - vol_0(ne_min-1) )/dvol_0(ne_min-1)
        fracti             = ( vol(i+1) - vol_0(ne_min-1) )/dvol_0(ne_min-1)
        rho(i)             = rhol(n) + half * ( fractim1 + fracti ) * ( drho(n) + rho6(n) ) &
&                          - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * rho6(n)
        dm(i)              = rho(i) * dvol(i)
        m(i+1)             = m(i) + dm(i)

      ELSE ! r(i+1) > r0(ne_min)

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i+1 < ne
!
!        :        |     :          :             :     |    :
!        :        |     :          :             :     |    :
!        :        |     :          :     ...     :     |    :
!        :        |     :          :             :     |    :
!        :        |     :          :             :     |    :
!
!     ne_min-1    i   ne_min   ne_min+1         ne-1 i+1   ne
!
!                 <--------dvolume_i_nem1------->
!                 <---------dmass_i_nem1-------->
!
!-----------------------------------------------------------------------
!  dmass_i_nem1 is initially the mass from i to ne_min
!  dvolume_i_nem1 is initially the volume from i to ne_min, calculated above
!-----------------------------------------------------------------------

        n                  = ne_min - 1 + 6
        d_vol_i_ne_min     = vol_0(ne_min) - vol(i)
        fractn             = half * d_vol_i_ne_min/dvol_0(ne_min-1)
        fractn2            = 1.d0 - fourthd * fractn
        rho_lft            = rhol(n) + drho(n) - fractn * ( drho(n) - fractn2 * rho6(n) )
        dvolume_i_nem1     = d_vol_i_ne_min
        dmass_i_nem1       = rho_lft * d_vol_i_ne_min
        DO nep = ne_min, i_out
          IF ( r(i+1) >= r0(nep+1) ) THEN
            dmass_i_nem1   = dmass_i_nem1 + dm0(nep)
            dvolume_i_nem1 = dvolume_i_nem1 + dvol_0(nep)
          ELSE ! r(i+1) < r0(nep+1)
            ne             = nep + 1
            nem1           = nep
            n              = nem1 + 6
            fractn         = half * ( vol(i+1) - vol_0(nem1) )/dvol_0(nem1)
            fractn2        = 1.d0 - fourthd * fractn
            rho_rgh        = rhol(n) + fractn * ( drho(n) + fractn2 * rho6(n) )
            dvolume_i_nem1 = dvolume_i_nem1 + vol(i+1) - vol_0(nem1)
            dmass_i_nem1   = dmass_i_nem1 + rho_rgh * ( vol(i+1) - vol_0(nem1) )
            EXIT
          END IF ! r(i+1) >= r0(nep+1)
        END DO ! nep = ne_min, i_out
        dm(i)              = dmass_i_nem1
        m(i+1)             = m(i) + dm(i)
        rho(i)             = dm(i)/dvol(i)

      END IF ! r(i+1) <= r0(ne_min)

    END DO Outer_i_loop2

!-----------------------------------------------------------------------
!                   ||||| Radial zones i_out |||||
!-----------------------------------------------------------------------

    dm(i_out)              = m0(i_out+1) - m(i_out)
    m(i_out+1)             = m0(i_out+1)
    rho(i_out)             = dm(i_out)/dvol(i_out)

!-----------------------------------------------------------------------
!              ||||| Radial zones i_out+1 to imax |||||
!-----------------------------------------------------------------------

    dm(i_out+1:imax)       = dm0(i_out+1:imax)
    m(i_out+2:imax+1)      = m0(i_out+2:imax+1)
    rho(i_out+1:imax)      = rho0(i_out+1:imax)

!-----------------------------------------------------------------------
!
!        \\\\\ PPM INTERPOLATE IN MASS THE TEMPERATURES, /////
!        \\\\\    ELECTRON FRACTIONS, AND VELOCITIES     /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Put pre-rezoned temperatures, electron fractions and velocities for
!   ray ij_ray, ik_ray into 1-D arrays
!-----------------------------------------------------------------------

    t0(:)                  = t_c (:,ij_ray,ik_ray)
    ye0(:)                 = ye_c(:,ij_ray,ik_ray)
    u0(:)                  = u_c (:,ij_ray,ik_ray)
    v0(:)                  = v_c (:,ij_ray,ik_ray)
    w0(:)                  = w_c (:,ij_ray,ik_ray)


!-----------------------------------------------------------------------
!  Put original mass variables into arrays, padding with 6 ghost zones
!
!  dxm(n) must be indexed as the shell above enxlosed mass xam(n)
!-----------------------------------------------------------------------

    xam0(nmin:nmax+1)      = m0(imin:imax+1)
    dxm0(nmin:nmax)        = dm0(imin:imax)

!-----------------------------------------------------------------------
!  Include six ghost zones on the left and right
!-----------------------------------------------------------------------

    DO n = 1, 6
      dxm0(nmin-n)         = dxm0(nmin+n-1)
      xam0 (nmin-n)        = xam0(nmin-n+1) - dxm0(nmin-n)
      dxm0(nmax+n)         = dxm0(nmax)
      xam0(nmax+n)         = xam0(nmax+n-1) + dxm0(nmax+n)
    END DO

!-----------------------------------------------------------------------
!  Compute parabolic coefficients in mass
!-----------------------------------------------------------------------

    CALL paraset( nx+12, param, dxm0, xam0, nmin-4, nmax+4, ngeom )

!-----------------------------------------------------------------------
!  Put original electron fractions, temperatures, and compositions
!   into arrays, padding with 6 ghost zones
!
!  t_pad(n), etc. must be indexed as the shells above enxlosed mass
!   xam(n)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Put original variables in padded arrays
!-----------------------------------------------------------------------

    t_pad (nmin:nmax)      = t0 (imin:imax)
    ye_pad(nmin:nmax)      = ye0(imin:imax)
    u_pad (nmin:nmax)      = u0 (imin:imax)
    v_pad (nmin:nmax)      = v0 (imin:imax)
    w_pad (nmin:nmax)      = w0 (imin:imax)

!-----------------------------------------------------------------------
!  include six ghost zones on the left
!-----------------------------------------------------------------------

    DO n = 1, 6
      t_pad (nmin-n)       = t_pad (nmin+n-1)
      ye_pad(nmin-n)       = ye_pad(nmin+n-1)
      u_pad (nmin-n)       = u_pad (nmin+n-1)
      v_pad (nmin-n)       = v_pad (nmin+n-1)
      w_pad (nmin-n)       = w_pad (nmin+n-1)
    END DO

!-----------------------------------------------------------------------
!  Generate parabolic interpolation coefficients
!-----------------------------------------------------------------------

    CALL parabola( nmin-3, nmax+3, ntot, param, t_pad , dt , t6 , tl , dum, &
&    0, 0, ngeom )
    CALL parabola( nmin-3, nmax+3, ntot, param, ye_pad, dye, ye6, yel, dum, &
&    0, 0, ngeom )
    CALL parabola( nmin-3, nmax+3, ntot, param, u_pad , du , u6 , ul , dum, &
&    0, 0, ngeom )
    CALL parabola( nmin-3, nmax+3, ntot, param, v_pad , dv , v6 , vl , dum, &
&    0, 0, ngeom )
    CALL parabola( nmin-3, nmax+3, ntot, param, w_pad , dw , w6 , wl , dum, &
&    0, 0, ngeom )

!-----------------------------------------------------------------------
!                   ||||| Radial zone 1 |||||
!-----------------------------------------------------------------------

    t (1)                  = t0 (1)
    ye(1)                  = ye0(1)
    u (1)                  = u0 (1)
    v (1)                  = v0 (1)
    w (1)                  = w0 (1)

!-----------------------------------------------------------------------
!                   ||||| Radial zone 2 |||||
!-----------------------------------------------------------------------

    ne_min                 = 2
    i                      = 2

!-----------------------------------------------------------------------
!  Find position of new zone edge relative to the nearby zone edges
!   in the original grid.
!
!  ne is the index of the original grid for which
!    r(ne-1) < r(i+1) < r(ne)
!-----------------------------------------------------------------------

    t_incl                 = zero
    ye_incl                = zero
    u_incl                 = zero
    v_incl                 = zero
    w_incl                 = zero
    DO nep = ne_min, i_out
      ne                   = nep + 1
      nem1                 = nep
      IF ( r0(nep+1) > r(i+1) ) EXIT
      t_incl               = t_incl  + t0(nem1)  * dm0(nem1)
      ye_incl              = ye_incl + ye0(nem1) * dm0(nem1)
      u_incl               = u_incl  + u0(nem1)  * dm0(nem1)
      v_incl               = v_incl  + v0(nem1)  * dm0(nem1)
      w_incl               = w_incl  + w0(nem1)  * dm0(nem1)
    END DO ! nep = ne_min, i_out

!-----------------------------------------------------------------------
!  Determine the mass and volume of material between nem1 and i
!
!        :             :     |    :
!        :             :     |    :
!        :     ...     :     |    :
!        :             :     |    :
!        :             :     |    :
!
!    ne_min=2         nem1  i+1   ne
!        i
!                      <---- n --->
!
!-----------------------------------------------------------------------

    n                      = nem1 + 6
    fractn                 = half * ( m(i+1) - m0(nem1) )/( m0(ne) - m0(nem1) )
    fractn2                = 1.d0 - fourthd * fractn

    t_rt                   = tl(n)   + fractn * ( dt(n)  + fractn2 * t6(n)  )
    t_incl                 = t_incl  + t_rt  * ( m(i+1) - m0(nem1) )

    ye_rt                  = yel(n)  + fractn * ( dye(n) + fractn2 * ye6(n) )
    ye_incl                = ye_incl + ye_rt * ( m(i+1) - m0(ne-1) )

    u_rt                   = ul(n)   + fractn * ( du(n)  + fractn2 * u6(n)  )
    u_incl                 = u_incl  + u_rt  * ( m(i+1) - m0(nem1) )

    v_rt                   = vl(n)   + fractn * ( dv(n)  + fractn2 * v6(n)  )
    v_incl                 = v_incl  + v_rt  * ( m(i+1) - m0(nem1) )

    w_rt                   = wl(n)   + fractn * ( dw(n)  + fractn2 * w6(n)  )
    w_incl                 = w_incl  + w_rt  * ( m(i+1) - m0(nem1) )

!-----------------------------------------------------------------------
!  Average
!-----------------------------------------------------------------------

    t (i)                  = t_incl /dm(i)
    ye(i)                  = ye_incl/dm(i)
    u (i)                  = u_incl /dm(i)
    v (i)                  = v_incl /dm(i)
    w (i)                  = w_incl /dm(i)

!-----------------------------------------------------------------------
!              ||||| Radial zones 3 < i < i_out - 1 |||||
!-----------------------------------------------------------------------

    Outer_i_loop3: DO i = 3, i_out - 1

      ne_min               = ne

      IF ( r(i+1) <= r0(ne_min) ) THEN

!-----------------------------------------------------------------------
!  Does i+1 lie between ne_min-1 and ne_min?
!
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    : 
!
!   ne_min-1    i         i+1  ne_min
!
!        <---------- n --------->
!
!  If so, average between i and i+1
!-----------------------------------------------------------------------

        n                = ne_min - 1 + 6
        fractim1         = ( m(i  ) - m0(ne_min-1) )/( m0(ne_min) - m0(ne_min-1) )
        fracti           = ( m(i+1) - m0(ne_min-1) )/( m0(ne_min) - m0(ne_min-1) )

        t(i)             = tl(n) + half * ( fractim1 + fracti ) * ( dt(n) + t6(n) ) &
&                        - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * t6(n)

        ye(i)            = yel(n) + half * ( fractim1 + fracti ) * ( dye(n) + ye6(n) ) &
&                        - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * ye6(n)

        u(i)             = ul(n) + half * ( fractim1 + fracti ) * ( du(n) + u6(n) ) &
&                        - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * u6(n)

        v(i)             = vl(n) + half * ( fractim1 + fracti ) * ( dv(n) + v6(n) ) &
&                        - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * v6(n)

        w(i)             = wl(n) + half * ( fractim1 + fracti ) * ( dw(n) + w6(n) ) &
&                        - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * w6(n)

      ELSE ! r(i+1) > r0(ne_min)

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i+1 < ne
!
!        :        |     :          :             :     |    :
!        :        |     :          :             :     |    :
!        :        |     :          :     ...     :     |    :
!        :        |     :          :             :     |    :
!        :        |     :          :             :     |    :
!
!     ne_min-1    i   ne_min   ne_min+1         ne-1 i+1   ne
!
!-----------------------------------------------------------------------

        t_incl             = zero
        ye_incl            = zero
        u_incl             = zero
        v_incl             = zero
        w_incl             = zero

        n                  = ne_min - 1 + 6
        dm_i_ne_min        = m0(ne_min) - m(i)
        fractn             = half * dm_i_ne_min/dm0(ne_min-1)
        fractn2            = 1.d0 - fourthd * fractn

        t_lt               = tl(n)  + dt(n)  - fractn * ( dt(n)  - fractn2 * t6(n)  )
        ye_lt              = yel(n) + dye(n) - fractn * ( dye(n) - fractn2 * ye6(n) )
        u_lt               = ul(n)  + du(n)  - fractn * ( du(n)  - fractn2 * u6(n)  )
        v_lt               = vl(n)  + dv(n)  - fractn * ( dv(n)  - fractn2 * v6(n)  )
        w_lt               = wl(n)  + dw(n)  - fractn * ( dw(n)  - fractn2 * w6(n)  )

        t_incl             = t_lt  * dm_i_ne_min
        ye_incl            = ye_lt * dm_i_ne_min
        u_incl             = u_lt  * dm_i_ne_min
        v_incl             = v_lt  * dm_i_ne_min
        w_incl             = w_lt  * dm_i_ne_min

        DO nep = ne_min, i_out
          IF ( r(i+1) >= r0(nep+1) ) THEN

            t_incl         = t_incl  + t0(nep)  * dm0(nep)
            ye_incl        = ye_incl + ye0(nep) * dm0(nep)
            u_incl         = u_incl  + u0(nep)  * dm0(nep)
            v_incl         = v_incl  + v0(nep)  * dm0(nep)
            w_incl         = w_incl  + w0(nep)  * dm0(nep)

          ELSE ! r(i+1) < r0(nep+1)

            ne             = nep + 1
            nem1           = nep
            n              = nem1 + 6
            fractn         = half * ( m(i+1) - m0(nem1) )/dm0(nem1)
            fractn2        = 1.d0 - fourthd * fractn

            t_rt           = tl(n)  + fractn * ( dt(n)  + fractn2 * t6(n)  )
            ye_rt          = yel(n) + fractn * ( dye(n) + fractn2 * ye6(n) )
            u_rt           = ul(n)  + fractn * ( du(n)  + fractn2 * u6(n)  )
            v_rt           = vl(n)  + fractn * ( dv(n)  + fractn2 * v6(n)  )
            w_rt           = wl(n)  + fractn * ( dw(n)  + fractn2 * w6(n)  )

            t_incl         = t_incl   + t_rt   * ( m(i+1) - m0(nem1) )
            ye_incl        = ye_incl  + ye_rt  * ( m(i+1) - m0(nem1) )
            u_incl         = u_incl   + u_rt   * ( m(i+1) - m0(nem1) )
            v_incl         = v_incl   + v_rt   * ( m(i+1) - m0(nem1) )
            w_incl         = w_incl   + w_rt   * ( m(i+1) - m0(nem1) )

            EXIT

          END IF ! r(i+1) >= r0(nep+1)
        END DO ! nep = ne_min, i_out

!-----------------------------------------------------------------------
!  Average
!-----------------------------------------------------------------------

        t (i)              = t_incl /dm(i)
        ye(i)              = ye_incl/dm(i)
        u (i)              = u_incl /dm(i)
        v (i)              = v_incl /dm(i)
        w (i)              = w_incl /dm(i)

      END IF ! r(i+1) <= r0(ne_min)

    END DO Outer_i_loop3

!-----------------------------------------------------------------------
!                   ||||| Radial zones i_out |||||
!-----------------------------------------------------------------------

    ne_min                 = ne

    IF ( r(i_out+1) <= r0(ne_min) ) THEN

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i+1 < ne
!
!        :         |                  |
!        :         |                  |
!        :         |                  |
!        :         |                  |
!        :         |                  |
!
!     ne_min-1   i_out             i_out+1
!                                  ne_min
!-----------------------------------------------------------------------

      n                    = ne_min - 1 + 6
      fractn               = half * dm(i_out)/dm0(ne_min-1)
      fractn2              = 1.d0 - fourthd * fractn

      t (i_out)            = tl(n)  + dt(n)  - fractn * ( dt(n)  - fractn2 * t6(n)  )
      ye(i_out)            = yel(n) + dye(n) - fractn * ( dye(n) - fractn2 * ye6(n) )
      u (i_out)            = ul(n)  + du(n)  - fractn * ( du(n)  - fractn2 * u6(n)  )
      v (i_out)            = vl(n)  + dv(n)  - fractn * ( dv(n)  - fractn2 * v6(n)  )
      w (i_out)            = wl(n)  + dw(n)  - fractn * ( dw(n)  - fractn2 * w6(n)  )

    ELSE ! r(i_out+1) > r0(ne_min)

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i+1 < ne
!
!        :         |          :          :             :     |
!        :         |          :          :             :     |
!        :         |          :          :     ...     :     |
!        :         |          :          :             :     |
!        :         |          :          :             :     |
!
!     ne_min-1   i_out     ne_min     ne_min+1       ne-1 i_out+1
!                                                            ne
!-----------------------------------------------------------------------

      t_incl               = zero
      ye_incl              = zero
      u_incl               = zero
      v_incl               = zero
      w_incl               = zero

      n                    = ne_min - 1 + 6
      dm_i_ne_min          = m0(ne_min) - m(i_out)
      fractn               = half * dm_i_ne_min/dm0(ne_min-1)
      fractn2              = 1.d0 - fourthd * fractn

      t_lt                 = tl(n)  + dt(n)  - fractn * ( dt(n)  - fractn2 * t6(n)  )
      ye_lt                = yel(n) + dye(n) - fractn * ( dye(n) - fractn2 * ye6(n) )
      u_lt                 = ul(n)  + du(n)  - fractn * ( du(n)  - fractn2 * u6(n)  )
      v_lt                 = vl(n)  + dv(n)  - fractn * ( dv(n)  - fractn2 * v6(n)  )
      w_lt                 = wl(n)  + dw(n)  - fractn * ( dw(n)  - fractn2 * w6(n)  )

      t_incl               = t_lt  * dm_i_ne_min
      ye_incl              = ye_lt * dm_i_ne_min
      u_incl               = u_lt  * dm_i_ne_min
      v_incl               = v_lt  * dm_i_ne_min
      w_incl               = w_lt  * dm_i_ne_min

      DO nep = ne_min, i_out
        IF ( r(i_out+1) >= r0(nep+1) ) THEN

          t_incl           = t_incl  + t0(nep)  * dm0(nep)
          ye_incl          = ye_incl + ye0(nep) * dm0(nep)
          u_incl           = u_incl  + u0(nep)  * dm0(nep)
          v_incl           = v_incl  + v0(nep)  * dm0(nep)
          w_incl           = w_incl  + w0(nep)  * dm0(nep)

        ELSE ! r(i_out+1) < r0(nep+1)

          EXIT

        END IF ! r(i_out+1) >= r0(nep+1)
      END DO ! nep = ne_min, i_out

!-----------------------------------------------------------------------
!  Average
!-----------------------------------------------------------------------

      t (i_out)            = t_incl /dm(i_out)
      ye(i_out)            = ye_incl/dm(i_out)
      u (i_out)            = u_incl /dm(i_out)
      v (i_out)            = v_incl /dm(i_out)
      w (i_out)            = w_incl /dm(i_out)

    END IF ! r(i+1) <= r0(ne_min)

!-----------------------------------------------------------------------
!              ||||| Radial zones i_out+1 to imax |||||
!-----------------------------------------------------------------------

    t (i_out+1:imax)       = t0 (i_out+1:imax)
    ye(i_out+1:imax)       = ye0(i_out+1:imax)
    u (i_out+1:imax)       = u0 (i_out+1:imax)
    v (i_out+1:imax)       = v0 (i_out+1:imax)
    w (i_out+1:imax)       = w0 (i_out+1:imax)

!-----------------------------------------------------------------------
!
!        \\\\\ PPM INTERPOLATE IN VOLUME THE ZERO-MOMENTS /////
!        \\\\\     OF THE NEUTRINO OCCUPATION NUMBER      /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Put pre-rezoned occupation number zero moments in 3-D arrays
!-----------------------------------------------------------------------

    psi0_0(:,:,:)          = psi0_c(:,:,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Load neutrino distribution function, padding with 6 ghost zones
!-----------------------------------------------------------------------

    psi0_re(nmin:nmax,:,:) = psi0_0(imin:imax,:,:)

!-----------------------------------------------------------------------
!  Load left (inner) ghosts
!-----------------------------------------------------------------------

    DO n = 1, 6
      nminn1               = MIN( nmin + n - 1, nmax )
      psi0_re(nmin-n,:,:)  = psi0_re(nminn1,:,:)
    END DO

!-----------------------------------------------------------------------
!  Generate interpolation functions for psi0.
!-----------------------------------------------------------------------

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)

        DO ii = 1,imax+12
          psi(ii)          = psi0_re(ii,k,nn)
        END DO !  ii = 1,imax+12

        CALL parabola( nmin-3, nmax+3, ntot, paraV, psi, dpsi, psi6, psil,  &
&        dm, 0, 0, ngeom )

        DO ii = 1,imax+12
          dpsi_a(ii,k,nn)  = dpsi(ii)
          psi6_a(ii,k,nn)  = psi6(ii)
          psil_a(ii,k,nn)  = psil(ii)
        END DO !  ii = 1,imax+12

      END DO ! k = 1,nnugp(nn)
    END DO ! nn = 1,nnu

!-----------------------------------------------------------------------
!                   ||||| Radial zone 1 |||||
!-----------------------------------------------------------------------

    psi0(1,:,:)            = psi0_0(1,:,:)

!-----------------------------------------------------------------------
!                   ||||| Radial zone 2 |||||
!-----------------------------------------------------------------------

    ne_min                 = 2
    i                      = 2

!-----------------------------------------------------------------------
!  Find position of new zone edge relative to the nearby zone edges
!   in the original grid.
!
!  ne is the index of the original grid for which
!    r(ne-1) < r(i+1) < r(ne)
!-----------------------------------------------------------------------

    psi_incl               = zero
    DO nep = ne_min, i_out
      ne                   = nep + 1
      nem1                 = nep
      IF ( r0(nep+1) > r(i+1) ) EXIT
      DO nn = 1,nnu
        IF ( nnugp(nn)  ==  0 ) CYCLE
        DO k = 1,nnugp(nn)
          psi_incl(k,nn)   = psi_incl(k,nn) + psi0_0(nem1,k,nn)  * dvol_0(nem1)
        END DO ! k = 1,nnugp(nn)
      END DO ! nn = 1,nnu
    END DO ! nep = ne_min, i_out

!-----------------------------------------------------------------------
!  Determine the mass and volume of material between nem1 and i
!
!        :             :     |    :
!        :             :     |    :
!        :     ...     :     |    :
!        :             :     |    :
!        :             :     |    :
!
!    ne_min=2         nem1  i+1   ne
!        i
!                      <---- n --->
!
!-----------------------------------------------------------------------

    n                      = nem1 + 6
    fractn                 = half * ( vol(i+1) - vol_0(nem1) )/dvol_0(nem1)
    fractn2                = 1.d0 - fourthd * fractn

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)
        psi_rt(k,nn)       = psil_a(n,k,nn) + fractn * ( dpsi_a(n,k,nn)  + fractn2 * psi6_a(n,k,nn) )
        psi_incl(k,nn)     = psi_incl(k,nn) + psi_rt(k,nn)  * ( vol(i+1) - vol_0(nem1) )
      END DO ! k = 1,nnugp(nn)
    END DO ! nn = 1,nnu

!-----------------------------------------------------------------------
!  Average
!-----------------------------------------------------------------------

    DO nn = 1,nnu
      IF ( nnugp(nn)  ==  0 ) CYCLE
      DO k = 1,nnugp(nn)
        psi0(i,k,nn)       = psi_incl(k,nn)/dvol(i)
      END DO ! k = 1,nnugp(nn)
    END DO ! nn = 1,nnu


!-----------------------------------------------------------------------
!              ||||| Radial zones 3 < i < i_out - 1 |||||
!-----------------------------------------------------------------------

    Outer_i_loop4: DO i = 3, i_out - 1

      ne_min               = ne

      IF ( r(i+1) <= r0(ne_min) ) THEN

!-----------------------------------------------------------------------
!  Does i+1 lie between ne_min-1 and ne_min?
!
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    :
!        :      |          |    : 
!
!   ne_min-1    i         i+1  ne_min
!
!        <---------- n --------->
!
!  If so, average between i and i+1
!-----------------------------------------------------------------------

        n                = ne_min - 1 + 6
        fractim1         = ( vol(i  ) - vol_0(ne_min-1) )/dvol_0(ne_min-1)
        fracti           = ( vol(i+1) - vol_0(ne_min-1) )/dvol_0(ne_min-1)

        DO nn = 1,nnu
          IF ( nnugp(nn)  ==  0 ) CYCLE
          DO k = 1,nnugp(nn)
            psi0(i,k,nn) = psil_a(n,k,nn) + half * ( fractim1 + fracti ) * ( dpsi_a(n,k,nn) + psi6_a(n,k,nn) ) &
&                        - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * psi6_a(n,k,nn)
          END DO ! k = 1,nnugp(nn)
        END DO ! nn = 1,nnu


      ELSE ! r(i+1) > r0(ne_min)

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i+1 < ne
!
!        :        |     :          :             :     |    :
!        :        |     :          :             :     |    :
!        :        |     :          :     ...     :     |    :
!        :        |     :          :             :     |    :
!        :        |     :          :             :     |    :
!
!     ne_min-1    i   ne_min   ne_min+1         ne-1 i+1   ne
!
!-----------------------------------------------------------------------

        psi_incl           = zero

        n                  = ne_min - 1 + 6
        d_vol_i_ne_min     = vol_0(ne_min) - vol(i)
        fractn             = half * d_vol_i_ne_min/dvol_0(ne_min-1)
        fractn2            = 1.d0 - fourthd * fractn

        DO nn = 1,nnu
          IF ( nnugp(nn)  ==  0 ) CYCLE
          DO k = 1,nnugp(nn)
            psi_lt(k,nn)   = psil_a(n,k,nn)  + dpsi_a(n,k,nn)  - fractn * ( dpsi_a(n,k,nn)  - fractn2 * psi6_a(n,k,nn)  )
            psi_incl(k,nn) = psi_lt(k,nn) * d_vol_i_ne_min
          END DO ! k = 1,nnugp(nn)
        END DO ! nn = 1,nnu

        DO nep = ne_min, i_out
          IF ( r(i+1) >= r0(nep+1) ) THEN

          DO nn = 1,nnu
            IF ( nnugp(nn)  ==  0 ) CYCLE
            DO k = 1,nnugp(nn)
              psi_incl(k,nn) = psi_incl(k,nn) + psi0_0(nep,k,nn) * dvol_0(nep)
            END DO ! k = 1,nnugp(nn)
          END DO ! nn = 1,nnu

          ELSE ! r(i+1) < r0(nep+1)

            ne             = nep + 1
            nem1           = nep
            n              = nem1 + 6
            fractn         = half * ( vol(i+1) - vol_0(nem1) )/dvol_0(nem1)
            fractn2        = 1.d0 - fourthd * fractn

            DO nn = 1,nnu
              IF ( nnugp(nn)  ==  0 ) CYCLE
              DO k = 1,nnugp(nn)
                psi_rt(k,nn)   = psil_a(n,k,nn)  + fractn * ( dpsi_a(n,k,nn)  + fractn2 * psi6_a(n,k,nn) )
                psi_incl(k,nn) = psi_incl(k,nn) + psi_rt(k,nn) * ( vol(i+1) - vol_0(nem1) )
              END DO ! k = 1,nnugp(nn)
            END DO ! nn = 1,nnu

            EXIT

          END IF ! r(i+1) >= r0(nep+1)
        END DO ! nep = ne_min, i_out

!-----------------------------------------------------------------------
!  Average
!-----------------------------------------------------------------------

        DO nn = 1,nnu
          IF ( nnugp(nn)  ==  0 ) CYCLE
          DO k = 1,nnugp(nn)
            psi0(i,k,nn)   = psi_incl(k,nn)/dvol(i)
          END DO ! k = 1,nnugp(nn)
        END DO ! nn = 1,nnu

      END IF ! r(i+1) <= r0(ne_min)

    END DO Outer_i_loop4

!-----------------------------------------------------------------------
!                   ||||| Radial zones i_out |||||
!-----------------------------------------------------------------------

    ne_min                 = ne

    IF ( r(i_out+1) <= r0(ne_min) ) THEN

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i+1 < ne
!
!        :         |                  |
!        :         |                  |
!        :         |                  |
!        :         |                  |
!        :         |                  |
!
!     ne_min-1   i_out             i_out+1
!                                  ne_min
!-----------------------------------------------------------------------

      n                    = ne_min - 1 + 6
      fractn               = half * dvol(i_out)/dvol_0(ne_min-1)
      fractn2              = 1.d0 - fourthd * fractn

      DO nn = 1,nnu
        IF ( nnugp(nn)  ==  0 ) CYCLE
        DO k = 1,nnugp(nn)
          psi0(i_out,k,nn) = psil_a(n,k,nn)  + dpsi_a(n,k,nn)             &
&                          - fractn * ( dpsi_a(n,k,nn)  - fractn2 * psi6_a(n,k,nn) )
        END DO ! k = 1,nnugp(nn)
      END DO ! nn = 1,nnu

    ELSE ! r(i_out+1) > r0(ne_min)

!-----------------------------------------------------------------------
!  FInd ne for which ne-1 < i+1 < ne
!
!        :         |          :          :             :     |
!        :         |          :          :             :     |
!        :         |          :          :     ...     :     |
!        :         |          :          :             :     |
!        :         |          :          :             :     |
!
!     ne_min-1   i_out     ne_min     ne_min+1       ne-1 i_out+1
!                                                            ne
!-----------------------------------------------------------------------

      psi_incl             = zero

      n                    = ne_min - 1 + 6
      d_vol_i_ne_min       = vol_0(ne_min) - vol(i_out)
      fractn               = half * d_vol_i_ne_min/dvol_0(ne_min-1)
      fractn2              = 1.d0 - fourthd * fractn

      DO nn = 1,nnu
        IF ( nnugp(nn)  ==  0 ) CYCLE
        DO k = 1,nnugp(nn)
          psi_lt(k,nn)     = psil_a(n,k,nn)  + dpsi_a(n,k,nn)             &
&                          - fractn * ( dpsi_a(n,k,nn)  - fractn2 * psi6_a(n,k,nn)  )
          psi_incl(k,nn)   = psi_lt(k,nn) * d_vol_i_ne_min
        END DO ! k = 1,nnugp(nn)
      END DO ! nn = 1,nnu

      DO nep = ne_min, i_out
        IF ( r(i_out+1) >= r0(nep+1) ) THEN

          DO nn = 1,nnu
            IF ( nnugp(nn)  ==  0 ) CYCLE
            DO k = 1,nnugp(nn)
              psi_incl(k,nn) = psi_incl(k,nn) + psi0_0(nep,k,nn) * dvol_0(nep)
            END DO ! k = 1,nnugp(nn)
          END DO ! nn = 1,nnu

        ELSE ! r(i_out+1) < r0(nep+1)

          EXIT

        END IF ! r(i_out+1) >= r0(nep+1)
      END DO ! nep = ne_min, i_out

!-----------------------------------------------------------------------
!  Average
!-----------------------------------------------------------------------

      DO nn = 1,nnu
        IF ( nnugp(nn)  ==  0 ) CYCLE
        DO k = 1,nnugp(nn)
          psi0(i_out,k,nn) = psi_incl(k,nn)/dvol(i_out)
        END DO ! k = 1,nnugp(nn)
      END DO ! nn = 1,nnu

    END IF ! r(i+1) <= r0(ne_min)

!-----------------------------------------------------------------------
!              ||||| Radial zones i_out+1 to imax |||||
!-----------------------------------------------------------------------

    psi0(i_out+1:imax,:,:) = psi0_0(i_out+1:imax,:,:)

!-----------------------------------------------------------------------
!
!        \\\\\ RESTORE REZONED VARIABLES TO MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

    rho_c (imin:imax,ij_ray,ik_ray)     = rho(imin:imax)
    t_c   (imin:imax,ij_ray,ik_ray)     = t  (imin:imax)
    ye_c  (imin:imax,ij_ray,ik_ray)     = ye (imin:imax)
    u_c   (imin:imax,ij_ray,ik_ray)     = u  (imin:imax)
    v_c   (imin:imax,ij_ray,ik_ray)     = v  (imin:imax)
    w_c   (imin:imax,ij_ray,ik_ray)     = w  (imin:imax)
    psi0_c(imin:imax,:,:,ij_ray,ik_ray) = psi0(imin:imax,:,:)

  END DO ! ik_ray = 1, ik_ray_dim
END DO ! ij_ray = 1, ij_ray_dim

x_ef(imin:imax+1)          = r(imin:imax+1)
dx_cf(imin:imax)           = x_ef(imin+1:imax+1) - x_ef(imin:imax)
x_cf(imin:imax)            = half * ( x_ef(imin+1:imax+1) + x_ef(imin:imax) )

!-----------------------------------------------------------------------
!  Record the regridding
!-----------------------------------------------------------------------

WRITE (nlog,123) time, DMAX1( time - t_bounce, zero ), ncycle

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

 5000 CONTINUE

DEALLOCATE (rhobar_rgrd, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhobar        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rhobar, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhobar_rgrd   '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (mbar, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'mbar          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dmbar, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmbar         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (m, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'm             '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dm, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (r, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r             '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (r0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r0            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (vol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vol           '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol          '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (vol_0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vol_0         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol_0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol_0        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (xaV0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xaV0          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dxV0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dxV0          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (paraV, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paraV         '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (m0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'm0            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dm0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0           '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (xam0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xam0          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dxm0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dxm0          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (param, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'param         '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (rho0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho0          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rho, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho           '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rho_pad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_pad       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rhol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhol          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rho6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho6          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (drho, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'drho          '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (t0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't0            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't             '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (t_pad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_pad         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (tl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tl            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (t6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't6            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dt            '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ye0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye0           '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye_pad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye_pad        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (yel, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yel           '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye6           '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye           '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (u0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u0            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (u, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u             '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (u_pad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u_pad         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ul, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ul            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (u6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u6            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (du, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'du            '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (v0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v0            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (v, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v             '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (v_pad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v_pad         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (vl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vl            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (v6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v6            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dv, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dv            '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (w0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w0            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (w, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w             '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (w_pad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w_pad         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (wl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wl            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (w6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w6            '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dw            '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (psi0_0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0_0        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psil, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dpsi, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi           '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psil_a, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil_a        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi6_a, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6_a        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dpsi_a, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi_a        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi_lt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi_lt        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi_rt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi_rt        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi_incl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi_incl      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dum, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dum           '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE regridder

      
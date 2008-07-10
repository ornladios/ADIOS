SUBROUTINE regrid_final_grid( imin, imax, nx, ij_ray_dim, ny, ik_ray_dim, &
& nz, nse_c, x_ef, dx_cf, x_cf, rhobar0, regrid, rho_regrid, grid_frac,   &
& int_pre_b, int_post_b, time, t_bounce )
!-----------------------------------------------------------------------
!
!    File:         regrid_final_grid
!    Module:       regrid_final_grid
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
!  nse_c          : nuclear statistical equilibrium flag
!  x_ef           : x grid zone face locations at cycle end [cm]
!  dx_cf          : x_ef(i+1) - x_ef(i) [cm]
!  x_cf           : x grid zone midpoint locations at cycle end [cm]
!  rhobar0        : angularly averaged density [g cm^{-3}]
!  regrid         : regrid option switch
!  rho_regrid     : regrid up to density rho_regrid or 2 zones behind shock whichever is larger [g cm^{-3}]
!  grid_frac      : fraction of a grid width grid is allowed to move per time step
!  int_pre_b      : number of cycles between successive regrids before bounce
!  int_post_b     : number of cycles between successive regrids after bounce
!  time           : the elapsed time since the initiation of the calculation
!  t_bounce       : the time of core bounce
!
!    Output arguments:
!  x_ef           : x grid zone face locations at cycle end [cm]
!  dx_cf          : x grid zone face locations at cycle end [cm]
!  x_cf           : x grid zone face locations at cycle end [cm]
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  cycle_module, edit_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero, half, one, third, frpith, frpi
USE physcnst_module, ONLY: pi, msolar

USE cycle_module, ONLY: ncycle, ncynu_trns, intnu_trns
USE edit_module, ONLY: nprint, nlog

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

INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: nse_c ! nuclear statistical equilibrium flag

INTEGER, INTENT(in)                  :: int_pre_b           ! number of cycles between successive regrids before bounce
INTEGER, INTENT(in)                  :: int_post_b          ! number of cycles between successive regrids after bounce

REAL(KIND=double), INTENT(in)                 :: rho_regrid ! regrid up to density rho_regrid or 5 zones behind shock, whichever is larger [g cm^{-3}]
REAL(KIND=double), INTENT(in)                 :: grid_frac  ! fraction of a grid width grid is allowed to move per time step
REAL(KIND=double), INTENT(in), DIMENSION(nx)  :: rhobar0    ! angularly averaged density [g cm^{-3}]
REAL(KIND=double), INTENT(in)                 :: time       ! the elapsed time since the initiation of the calculation
REAL(KIND=double), INTENT(in)                 :: t_bounce   ! the time of core bounce

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx) :: dx_cf      ! radial grid widths [cm]
REAL(KIND=double), INTENT(out), DIMENSION(nx) :: x_cf       ! radial grid zone midpoints

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx+1) :: x_ef     ! radial grid edges [cm]

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=14)                      :: var_name

LOGICAL                                 :: first = .true.
LOGICAL                                 :: l_repeat = .false.
LOGICAL                                 :: l_bounds

INTEGER                                 :: istat            ! allocation status

INTEGER                                 :: it               ! iteration index
INTEGER                                 :: it_count         ! iteration counter
INTEGER, PARAMETER                      :: it_max = 50      ! maximum number of iteration attempts
INTEGER                                 :: i_nonconv        ! number of consecutive failure to converge in density

INTEGER, PARAMETER                      :: ngeom = 0        ! geometry glag

INTEGER                                 :: i                ! radial zone index
INTEGER                                 :: ii               ! radial zone index
INTEGER                                 :: ij_ray           ! j-index of a radial ray
INTEGER                                 :: ik_ray           ! k-index of a radial ray
INTEGER                                 :: n                ! padded index array
INTEGER                                 :: nmin             ! minimum padded radial array index
INTEGER                                 :: nmax             ! maximum padded radial array index
INTEGER                                 :: ntot             ! nmax + 6

INTEGER                                 :: i_in             ! innermost radial zone index of region to be rezoned
INTEGER                                 :: i_in_d           ! innermost radial zone index of region to be rezoned by density cutoff criteria
INTEGER, PARAMETER                      :: i_in_max = 30    ! maximum number of zones interior to region to be regridded
INTEGER                                 :: i_out            ! outermost radial zone index of region to be rezoned

INTEGER                                 :: i_10_hi          ! exponent of the upper density scale height in the region to be re-gridded
INTEGER                                 :: i_10_lo          ! exponent of the lower density  scale height in the region to be re-gridded
INTEGER                                 :: n_scale          ! number of density scale heights in the region to be re-gridded
INTEGER, DIMENSION(100)                 :: i_bounds         ! indices of the density interval boundaries

INTEGER, DIMENSION(nx)                  :: i_close          ! index of r0 closest to a given r

INTEGER                                 :: ne               ! X_old(ne) > X_new(i)
INTEGER                                 :: nem1             ! ne - 1
INTEGER                                 :: nep              ! do index
INTEGER                                 :: ne_min           ! minimum index in the new grid to search for the next zone edge

INTEGER                                 :: jjshockmx        ! maximum radial zone index of shock for a given ray
INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim)  :: jjshockmn        ! minimum radial zone index of shock for a given ray
INTEGER, DIMENSION(ny,nz)                  :: jjshockmn_a      ! minimum radial zone index of shock for a given ray
INTEGER                                 :: ishockmin        ! minimum radial zone index
INTEGER                                 :: ishockmin_prev = 1000 ! minimum radial zone index in previous regrid
INTEGER, DIMENSION(ij_ray_dim,ik_ray_dim)  :: inse_mn          ! outermost radial zone index for which nes = 1
INTEGER, DIMENSION(ny,nz)               :: inse_mn_a        ! outermost radial zone index for which nes = 1
INTEGER                                 :: inse_min         ! minimum radial zone index

REAL(KIND=double), PARAMETER                     :: tol = 1.d-6     ! iteration tolerence
REAL(KIND=double), PARAMETER                     :: pqmin = 2.d-1   ! minimum ratio of pseudoviscosity to pressure to locate shock

REAL(KIND=double)                                :: thpifrth        ! 3/4pi

REAL(KIND=double), PARAMETER                     :: rhobar_max = 1.d+14 ! if rhobar0(1) > rhobar_max re-grid only if time - t_bounce > t_b_min
REAL(KIND=double), PARAMETER                     :: rho_min_max = 1.d+15 ! re-grid only if rho_min < rho_min_max
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
REAL(KIND=double), DIMENSION(2000)               :: r0_prev         ! zone radii of the model after previous regrid [cm]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dr0             ! change in the model during preceding time step [cm]
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: r               ! radial grid achieved when maximum change criteria have benn appled
REAL(KIND=double), DIMENSION(2000)               :: r_new_grid      ! rezoned radial grid [cm]
REAL(KIND=double)                                :: delta_r         ! zone width of the region above 1d14 g cm^{-3} when regridded
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
REAL(KIND=double)                                :: dmass0          ! mass of one or multiple old shells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: m0              ! unshifted array of old enclosed masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dm0             ! unshifted array of old shell masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: xam0            ! padded array of old enclosed masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dxm0            ! padded array of old zone masses
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)   :: param           ! parabolic coefficients for masses

REAL(KIND=double)                                :: rho_test        ! dmass0/dvolume0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rho_pad         ! padded density array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rhol            ! density at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: rho6            ! density parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: drho            ! density slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)     :: dum             ! dummy array

REAL(KIND=double)                                :: rho_lft         ! left part of average density
REAL(KIND=double)                                :: rho_rgh         ! right part of average density

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Cannot find i_out in subroutine regrid_final_grid; rhobar0(1)=',es11.3, &
&  ' rhobar0(imax)=',es11.3,' rho_min=',es11.3,' run stopped')
  103 FORMAT (' nse(',i4,',',i4,',',i4,')=',i4,' cannot regrid non-NSE regions')
  105 FORMAT (' i_in_d=',i4,' > i_in_max=',i4,' i_in=',i4)
  111 FORMAT (' Iteration will not converge for rho_test in subroutine regrid_final_grid for determining zone ',i4,/ &
&      ' dvolume_test=',es11.3,' dvol_0(nem1)=',es11.3,' rho_test=',es11.3, &
&      ' rhobar_rgrd(i)=',es11.3)
  113 FORMAT (' Iteration will not converge for rho_test in subroutine regrid_final_grid for determining zone ',i4,/ &
&      ' for the case in which i and i-1 are in the same original zone',/ &
&      ' dvolume_test=',es11.3,' dvol_0(nem1)=',es11.3,' rho_test=',es11.3, &
&      ' rhobar_rgrd(i)=',es11.3,' run stopped')
  115 FORMAT (' Iteration will not converge for rho_test in subroutine regrid_final_grid for determining zone ',i4,/ &
&      ' for the case in which i and i-1 are in the same original zone',/ &
&      ' dvolume_test=',es11.3,' dvolume_i_nem1=',es11.3,' rho_test=',es11.3, &
&      ' rhobar_rgrd(i)=',es11.3)
  121 FORMAT (' Grid adjusted; i_out=',i4)
  123 FORMAT (' Regridding is being attempted in regrid_final_grid from rho_max=',es11.3, &
& ' to rho_min=',es11.3,' ishockmin=',i4,' inse_min=',i4,' i_in=',i4,' i_out=',i4,' rhobar0_shk=',es11.3, &
& ' rho_regrid=',es11.3) 
 1001 FORMAT (' Allocation problem for array ',a10,' in regrid')
 2001 FORMAT (' Deallocation problem for array ',a10,' in regrid')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if regrid = 'no'
!-----------------------------------------------------------------------

IF ( regrid == 'no' ) RETURN

!-----------------------------------------------------------------------
!  Return if ncycle <= 1
!-----------------------------------------------------------------------

IF ( ncycle <= 1 ) THEN
  first                 = .false.
  RETURN
END IF ! ncycle <= 1

!-----------------------------------------------------------------------
!
!       \\\\\ LET OLD GRID APPROACH NEW IF L_REPEAT = TRUE /////
!
!-----------------------------------------------------------------------

IF ( l_repeat ) THEN

  WRITE (nlog,121) i_out

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

  ALLOCATE (r(nx+1), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'r             '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (r0(nx+1), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'r0            '; WRITE (nlog,1001) var_name; END IF
  ALLOCATE (dr0(nx+1), STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'dr0           '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!   ||||| Move the rezoned grid with the fluid if m_grid is on |||||
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Move the rezoned grid with the moving grid
!  dr0 : change in the grid position due to the moving grid option
!-----------------------------------------------------------------------

  r0(:)                  = x_ef(:)
  dr0(1:imax+1)          = r0(1:imax+1) - r0_prev(1:imax+1)
  r(1:i_in)              = x_ef(1:i_in)
  DO i = 3, i_out + 1
    i_close(i)           = MINLOC( DABS( r0 - r_new_grid(i) ), DIM = 1 )
  END DO ! i = 3, i_out + 1

  DO i = 3, i_out + 1
    n                    = i_close(i) 
    r(i)                 = r_new_grid(i) + rinterp( r0(n-1), r0(n), r0(n+1), &
&    dr0(n-1), dr0(n), dr0(n+1), r_new_grid(i) )
  END DO ! i = 3, i_out + 1
  
  r(i_out+2:imax+1)      = x_ef(i_out+2:imax+1)

!-----------------------------------------------------------------------
!  Save rezoned radial grid
!-----------------------------------------------------------------------

  r_new_grid(1:imax+1)   = r(1:imax+1)

!-----------------------------------------------------------------------
!     ||||| Now move the actual grid towards the rezoned grid |||||
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Now limit change in r(i) to grid_frac * MIN( dr(i), dr(i-1) )
!-----------------------------------------------------------------------

  l_repeat               = .false.
  DO i = 3, i_out
    IF ( DABS( r(i) - r0(i) ) > grid_frac * DMIN1( r0(i+1) - r0(i), r0(i) - r0(i-1) ) ) THEN
      r(i)               = r0(i) + SIGN( grid_frac * DMIN1( r0(i+1) - r0(i), r0(i) - r0(i-1) ), r(i) - r0(i) )
      l_repeat           = .true.
    END IF ! r(i) - r0(i) > grid_frac * DMAX1( dr(i), dr(i-1) )
  END DO ! 2, i_out

!-----------------------------------------------------------------------
!  Save rezoned radial grid after maximum change criteria imposed
!-----------------------------------------------------------------------

  r0_prev(1:imax+1)      = r(1:imax+1)

!-----------------------------------------------------------------------
!  Return new grid
!-----------------------------------------------------------------------

  x_ef(imin:imax+1)      = r(imin:imax+1)
  dx_cf(imin:imax)       = x_ef(imin+1:imax+1) - x_ef(imin:imax)
  x_cf(imin:imax)        = half * ( x_ef(imin+1:imax+1) + x_ef(imin:imax) )

!-----------------------------------------------------------------------
!  Activate neutrino transport if high density zones have been adjusted
!-----------------------------------------------------------------------

  IF ( i_in /= i_in_d ) ncynu_trns = intnu_trns

!-----------------------------------------------------------------------
!  Deallocate arrays and return
!-----------------------------------------------------------------------

  DEALLOCATE (r, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'r             '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (r0, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'r0            '; WRITE (nlog,2001) var_name; END IF
  DEALLOCATE (dr0, STAT = istat)
    IF ( istat /= 0 ) THEN; var_name = 'dr0           '; WRITE (nlog,2001) var_name; END IF

END IF ! l_repeat

!-----------------------------------------------------------------------
!
!               \\\\\ DETERMINE WHETHER TO REGRID /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  If rhobar0(1) < rhobar_max, return if cycle number is not a multiple
!   of int_pre_b and l_repeat is false
!-----------------------------------------------------------------------

IF ( .not. first ) THEN
  IF ( rhobar0(1) < rhobar_max ) THEN
    IF ( MOD( ncycle, int_pre_b  ) /= 0 ) RETURN
  ELSE
    IF ( MOD( ncycle, int_post_b  ) /= 0 ) RETURN
  END IF ! rhobar0 > rhobar_max
END IF ! first

!-----------------------------------------------------------------------
!
!          \\\\\ DETERMINE DENSITY FROM WHICH TO REGRID /////
!
!-----------------------------------------------------------------------

i_in_d                   = MAXLOC( rhobar0, DIM = 1, MASK = rhobar0 < 1.d+14 )
i_in_d                   = MAX( i_in_d, 2 )
i_in                     = i_in_d
IF ( i_in > i_in_max ) THEN
  i_in                   = i_in - 1
  WRITE (nlog,105) i_in_d, i_in_max, i_in
END IF ! i_in > i_in_max

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

jjshockmn_a(:,:)         = jjshockmn(:,:)
inse_mn_a(:,:)           = inse_mn(:,:)
ishockmin                = MINVAL( jjshockmn_a(:,:) ) - 2
ishockmin                = MIN( ishockmin, ishockmin_prev + 1 )
ishockmin_prev           = ishockmin
inse_min                 = MINVAL( inse_mn_a(:,:) )
i_out                    = MIN( ishockmin, inse_min )
rhobar0_shk              = rhobar0(i_out)
rho_min                  = DMAX1( rhobar0_shk, rho_regrid )

!-----------------------------------------------------------------------
!  Return if rho_min > rho_min_max
!-----------------------------------------------------------------------

IF ( rho_min > rho_min_max ) RETURN

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
    i_out                = i
    EXIT
  END IF ! rhobar0(i) >= rho_min  .and.  rhobar0(i+1) <= rho_min
END DO

IF ( i_out == imax ) THEN
  WRITE (nprint,101) rhobar0(1), rhobar0(imax), rho_min
  WRITE (nlog,101) rhobar0(1), rhobar0(imax), rho_min
  STOP
END IF ! i_out == imax

DO ij_ray = 1, ij_ray_dim
  DO ik_ray = 1, ik_ray_dim
    IF ( nse_c(i_out,ij_ray,ik_ray) == 0 ) THEN
      RETURN
    END IF ! nse_c(i_out,ij_ray,ik_ray) == 0
  END DO ! ik_ray = 1, ik_ray_dim
END DO ! ij_ray = 1, ij_ray_dim

nmin                     = 7
nmax                     = imax + 6
ntot                     = nmax + 6

!-----------------------------------------------------------------------
!  Return if i_out <= i_in + 4
!-----------------------------------------------------------------------

IF ( i_out <= i_in + 4 ) RETURN

!-----------------------------------------------------------------------
!
!         \\\\\ DETERMINE R(I) GIVING CONSTANT DRHO/RHO /////
!
!-----------------------------------------------------------------------

WRITE (nlog,123) rhobar0(i_in_d), rho_min, ishockmin, inse_min, i_in, i_out, &
& rhobar0_shk, rho_regrid

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
  IF ( istat /= 0 ) THEN; var_name = 'r             '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (vol(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vol           '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol          '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (r0(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r0            '; WRITE (nlog,1001) var_name; END IF
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

ALLOCATE (rho_pad(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_pad       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rhol(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhol          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rho6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho6          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (drho(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'drho          '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dum(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dum           '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

i_nonconv                = 0

fourthd                  = 4.d0/3.d0
thpifrth                 = 1.d0/frpith

rhobar_rgrd              = zero
rhobar                   = zero

mbar                     = zero
dmbar                    = zero
m                        = zero
dm                       = zero
r                        = zero
vol                      = zero
dvol                     = zero

r0                       = zero
vol_0                    = zero
dvol_0                   = zero
xaV0                     = zero
dxV0                     = zero
paraV                    = zero

m0                       = zero
dm0                      = zero
xam0                     = zero
dxm0                     = zero
param                    = zero

rho_pad                  = zero
rhol                     = zero
rho6                     = zero
drho                     = zero

dum                      = zero

r0(:)                    = x_ef(:)

!-----------------------------------------------------------------------
!  Calculate original zone volumes
!-----------------------------------------------------------------------

vol_0(1)                 = zero
vol_0(2:imax+1)          = frpith * x_ef(2:imax+1)**3
dvol_0(1:imax)           = vol_0(2:imax+1) - vol_0(1:imax)

!-----------------------------------------------------------------------
!  Calculate original mass shells and enclosed masses
!-----------------------------------------------------------------------

dm0(1:imax)              = rhobar0(1:imax) * dvol_0(1:imax)

m0(1)                    = zero
DO i = 1,imax
  m0(i+1)                = m0(i) + dm0(i)
END DO ! 1,imax

!-----------------------------------------------------------------------
!  Now calculate the mean densities in the regridded region
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  rho_ratio is the ratio between adjacent zones to be achieved in the
!   regridded region
!
!              i_out - 1
!     rho_ratio          = rhobar0(i_out)/rhobar0(i_in_d-1)
!-----------------------------------------------------------------------

rhobar_rgrd(1:i_in-2)    = rhobar0(1:i_in-2)
rhobar_rgrd(i_in-1)      = rhobar0(i_in_d-1)

rho_ratio                = ( rhobar0(i_out)/rhobar0(i_in_d-1) )**( 1.d0/DBLE(i_out-i_in-1+1) )
DO i = i_in, i_out
  rhobar_rgrd(i)         =  rhobar_rgrd(i-1) * rho_ratio
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

xaV0(nmin:nmax+1)        = vol_0(imin:imax+1)
dxV0(nmin:nmax)          = dvol_0(imin:imax)

!-----------------------------------------------------------------------
!  Fill ghost zones
!-----------------------------------------------------------------------

V_ratio_b                = dxV0(nmax)/dxV0(nmax-1)
DO n = 1, 6
  dxV0(nmin-n)           = dxV0(nmin+n-1)
  xaV0(nmin-n)           = xaV0(nmin-n+1) - dxV0(nmin-n)
  dxV0(nmax+n)           = dxV0(nmax) * V_ratio_b**(n)
  xaV0(nmax+n)           = xaV0(nmax+n-1) + dxV0(nmax+n)
END DO

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

rho_pad(nmin:nmax)       = rhobar0(imin:imax)

!-----------------------------------------------------------------------
!  Fill ghost zones
!-----------------------------------------------------------------------

rho_ratio_b              = rho_pad(nmax)/rho_pad(nmax-1)
DO n = 1, 6
  rho_pad(nmin-n)        = rho_pad(nmin+n-1)
  rho_pad(nmax+n)        = rho_pad(nmax) * rho_ratio_b**(n)
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
!                 ||||| Radial zone 1 to i_in-1 |||||
!-----------------------------------------------------------------------

rhobar(1:i_in-1)         = rhobar0(1:i_in-1)
mbar(1)                  = zero
dmbar(1:i_in-1)          = dm0(1:i_in-1)
mbar(2:i_in)             = mbar(1:i_in-1) + dmbar(1:i_in-1)

r(1)                     = zero
r(2:i_in)                = r0(2:i_in)
IF ( i_in /= i_in_d ) THEN
  delta_r                = ( r0(i_in_d) - r0(2) )/DBLE( i_in - 2 )
  DO i = 3, i_in
    r(i)                 = r(i-1) + delta_r
  END DO ! i = 3, i_in
END IF ! i_in /= i_in_d

vol(1)                   = zero
vol(2:i_in)              = frpith * r(2:i_in)**3
dvol(1:i_in-1)           = vol(2:i_in) - vol(1:i_in-1)

!-----------------------------------------------------------------------
!                     ||||| Radial zone i_in |||||
!-----------------------------------------------------------------------

ne_min                   = i_in
i                        = i_in

!-----------------------------------------------------------------------
!  Find position of new zone edge relative to the nearby zone edges
!   in the original grid.
!
!  ne is the index of the original grid for which the mean density
!   between i and ne falls below rhobar_rgrd(i)
!-----------------------------------------------------------------------

dmass_2_nem1             = zero
dvolume_2_nem1           = zero
DO nep = ne_min, i_out
  dmass_2_nem1           = dmass_2_nem1 + dm0(nep)
  dvolume_2_nem1         = dvolume_2_nem1 + dvol_0(nep)
  rho_test               = dmass_2_nem1/dvolume_2_nem1
  IF ( rho_test < rhobar_rgrd(i) ) THEN
    ne                   = nep + 1
    nem1                 = nep
    dmass_2_nem1         = dmass_2_nem1 - dm0(nep)
    dvolume_2_nem1       = dvolume_2_nem1 - dvol_0(nep)
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
!   ne_min=i_in       nem1  i+1   ne
!        i
!                      <---- n --->
!
!-----------------------------------------------------------------------

n                        = nem1 + 6
dvolume_min              = zero
dvolume_max              = dvol_0(nem1)
DO it = 1, it_max

  it_count               = it
  dvolume_test           = half * ( dvolume_min + dvolume_max )
  fractn                 = half * dvolume_test/dvol_0(nem1)
  fractn2                = 1.d0 - fourthd * fractn
  rho_rgh                = rhol(n) + fractn * ( drho(n) + fractn2 * rho6(n) )
  rho_test               = ( dmass_2_nem1 + rho_rgh * dvolume_test )/( dvolume_2_nem1 + dvolume_test )
    
  IF ( DABS( rho_test - rhobar_rgrd(i) )/rhobar_rgrd(i) < tol ) EXIT

  IF ( rho_test >= rhobar_rgrd(i) ) THEN
    dvolume_min          = dvolume_test
  ELSE ! rho_test < rhobar_rgrd(i)
    dvolume_max          = dvolume_test
  END IF ! rho_test >= rhobar_rgrd(i)

END DO ! it = 1, it_max

IF ( it_count == it_max ) THEN
  WRITE (nprint,111) i, dvolume_test, dvol_0(nem1), rho_test, rhobar_rgrd(i)
  WRITE (nlog,111) i, dvolume_test, dvol_0(nem1), rho_test, rhobar_rgrd(i)
END IF ! it == it_max
  
!-----------------------------------------------------------------------
!  Success! Find new shell mass, enclosed mass, etc. of the new mass
!   shell
!-----------------------------------------------------------------------
  
rhobar(i)                = rho_test
dmbar(i)                 = dmass_2_nem1 + rho_rgh * dvolume_test
mbar(i+1)                = mbar(i) + dmbar(i)
r(i+1)                   = ( r(i)**3 + thpifrth * ( dvolume_2_nem1 + dvolume_test ) )**third
vol(i+1)                 = frpith * r(i+1)**3
dvol(i)                  = vol(i+1) - vol(i)

!-----------------------------------------------------------------------
!              ||||| Radial zones 3 < i < i_out - 1 |||||
!-----------------------------------------------------------------------

Outer_i_loop1: DO i = i_in + 1, i_out - 1

!-----------------------------------------------------------------------
!  ne_min -1 and ne_min are the original grid edges that are on the
!   immediate left and right, respectively, of new grid edge i
!
!  After each loop over i, ne lies just to the right of i
!-----------------------------------------------------------------------
  DO nep = 1, i_out
    IF ( r0(nep) >= r(i) ) THEN
      ne                 = nep
      nem1               = nep - 1
      EXIT
    END IF !  r0(nep) >= r(ii)
  END DO ! nep = 1, i_out

  ne_min                 = ne
  dmass0                 = zero
  dvolume0               = zero

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

  n                      = ne_min - 1 + 6
  d_vol_i_ne_min         = vol_0(ne_min) - vol(i)
  fractn                 = half * d_vol_i_ne_min/dvol_0(ne_min-1)
  fractn2                = 1.d0 - fourthd * fractn
  rho_lft                = rhol(n) + drho(n) - fractn * ( drho(n) - fractn2 * rho6(n) )

!-----------------------------------------------------------------------
!  If rho_lft <= rhobar_rgrd(i), i+1 lies between ne_min-1 and ne_min
!   to the right of i
!-----------------------------------------------------------------------

  IF ( rho_lft <= rhobar_rgrd(i) ) THEN

    it_count             = 0
    dvolume_min          = zero
    dvolume_max          = d_vol_i_ne_min
    DO it = 1, it_max

      it_count           = it
      dvolume_test       = half * ( dvolume_min + dvolume_max )
      fractim1           = ( vol(i)                - vol_0(ne_min-1) )/dvol_0(ne_min-1)
      fracti             = ( vol(i) + dvolume_test - vol_0(ne_min-1) )/dvol_0(ne_min-1)
      rho_test           = rhol(n) + half * ( fractim1 + fracti ) * ( drho(n) + rho6(n) ) &
&                        - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * rho6(n)
  
      IF ( DABS( rho_test - rhobar_rgrd(i) )/rhobar_rgrd(i) < tol ) EXIT

      IF ( rho_test >= rhobar_rgrd(i) ) THEN
        dvolume_min      = dvolume_test
      ELSE ! rho_test < rhobar_rgrd(i)
        dvolume_max      = dvolume_test
      END IF ! rho_test >= rhobar_rgrd(i)

    END DO ! it = 1, it_max

    IF ( it_count == it_max ) THEN
      i_nonconv          = i_nonconv + 1

!-----------------------------------------------------------------------
!  If rho_lft <= DSQRT( rhobar_rgrd(i) * rhobar_rgrd(i-1) ), temporarly
!   locate i+1 at ne_min
!-----------------------------------------------------------------------

      IF ( rho_lft <= DSQRT( rhobar_rgrd(i) * rhobar_rgrd(i-1) ) ) THEN
        rho_test         = rho_lft
        dvolume_test     = vol_0(ne_min) - vol(i)
      ELSE ! rho_lft > DSQRT( rhobar_rgrd(i) * rhobar_rgrd(i-1) )
        WRITE (nprint,113) i, dvolume_test, d_vol_i_ne_min, rho_test, rhobar_rgrd(i)
        WRITE (nlog,113) i, dvolume_test, d_vol_i_ne_min, rho_test, rhobar_rgrd(i)
        STOP
      END IF ! rho_lft <= DSQRT( rhobar_rgrd(i) * rhobar_rgrd(i-1) )
    END IF ! it_count == it_max
  
!-----------------------------------------------------------------------
!  Success! Find new shell mass, enclosed mass, etc. of the new mass
!   shell
!-----------------------------------------------------------------------
  
    rhobar(i)            = rho_test
    dmbar(i)             = rho_test * dvolume_test
    mbar(i+1)            = mbar(i) + dmbar(i)
    r(i+1)               = ( r(i)**3 + thpifrth * dvolume_test )**third
    vol(i+1)             = frpith * r(i+1)**3
    dvol(i)              = vol(i+1) - vol(i)

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

    dvolume_i_nem1       = d_vol_i_ne_min
    dmass_i_nem1         = rho_lft * d_vol_i_ne_min

    DO nep = ne_min, i_out

      dmass_i_nem1       = dmass_i_nem1 + dm0(nep)
      dvolume_i_nem1     = dvolume_i_nem1 + dvol_0(nep)
      rho_test           = dmass_i_nem1/dvolume_i_nem1

      IF ( rho_test < rhobar_rgrd(i) ) THEN
        ne               = nep + 1
        nem1             = nep
        dmass_i_nem1     = dmass_i_nem1 - dm0(nep)
        dvolume_i_nem1   = dvolume_i_nem1 - dvol_0(nep)
        EXIT
      END IF

    END DO ! nep = ne_min, i_out
  
!-----------------------------------------------------------------------
!  To determine r(i+1), iterate on the volume between r(i+1) and
!   r0(ne-1) until rho_test = rhobar_rgrd(i)
!-----------------------------------------------------------------------

    it_count             = 0
    n                    = nem1 + 6
    dvolume_min          = zero
    dvolume_max          = dvol_0(nem1)
    DO it = 1, it_max

      it_count           = it
      dvolume_test       = half * ( dvolume_min + dvolume_max )
      fractn             = half * dvolume_test/dvol_0(nem1)
      fractn2            = 1.d0 - fourthd * fractn
      rho_rgh            = rhol(n) + fractn * ( drho(n) + fractn2 * rho6(n) )
      rho_test           = ( dmass_i_nem1 + rho_rgh * dvolume_test )/( dvolume_i_nem1 + dvolume_test )
  
      IF ( DABS( rho_test - rhobar_rgrd(i) )/rhobar_rgrd(i) < tol ) EXIT

      IF ( rho_test >= rhobar_rgrd(i) ) THEN
        dvolume_min      = dvolume_test
      ELSE ! rho_test < rhobar_rgrd(i)
        dvolume_max      = dvolume_test
      END IF ! rho_test >= rhobar_rgrd(i)

    END DO ! it = 1, it_max

    IF ( it_count == it_max ) THEN
      i_nonconv          = i_nonconv + 1
      WRITE (nlog,115) i, dvolume_test, dvolume_i_nem1, rho_test, rhobar_rgrd(i)
    END IF ! it == it_max
  
!-----------------------------------------------------------------------
!  Success! Find new shell mass, enclosed mass, etc. of the new mass
!   shell
!-----------------------------------------------------------------------
  
    rhobar(i)            = rho_test
    dmbar(i)             = dmass_i_nem1 + rho_rgh * dvolume_test
    mbar(i+1)            = mbar(i) + dmbar(i)
    r(i+1)               = ( r(i)**3 + thpifrth * ( dvolume_i_nem1 + dvolume_test ) )**third
    vol(i+1)             = frpith * r(i+1)**3
    dvol(i)              = vol(i+1) - vol(i)

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
    r_ratio              = ( r(i+1)/r(i-i_nonconv) )**( 1.d0/DBLE(i_nonconv+1) )
    DO ii = 1, i_nonconv
      r(i-i_nonconv+ii)  = r(i-i_nonconv+ii-1) * r_ratio
      vol(i-i_nonconv+ii)= frpith * r(i-i_nonconv+ii)**3
      dvol(i-i_nonconv+ii-1) = vol(i-i_nonconv+ii) - vol(i-i_nonconv+ii-1)
    END DO ! ii = 1, i_nonconv
    dvol(i)              = vol(i+1) - vol(i)

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
        IF ( r0(nep) > r(ii) ) THEN
          ne             = nep
          nem1           = nep - 1
          EXIT
        END IF !  r0(nep) >= r(ii)
      END DO ! nep = 1, i_out
          
      ne_min             = ne
      dmass0             = zero
      dvolume0           = zero

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

        n                = ne_min - 1 + 6
        fractim1         = ( vol(ii  ) - vol_0(ne_min-1) )/dvol_0(ne_min-1)
        fracti           = ( vol(ii+1) - vol_0(ne_min-1) )/dvol_0(ne_min-1)
        rhobar(ii)       = rhol(n) + half * ( fractim1 + fracti ) * ( drho(n) + rho6(n) ) &
&                        - third * ( fractim1 * ( fractim1 + fracti ) + fracti * fracti ) * rho6(n)
        dmbar(ii)        = rhobar(ii) * dvol(ii)
        mbar(ii+1)       = mbar(ii) + dmbar(ii)

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

dmbar(i_out)             = m0(i_out+1) - mbar(i_out)
vol(i_out+1)             = vol_0(i_out+1)
dvol(i_out)              = vol(i_out+1) - vol(i_out)
rhobar(i_out)            = dmbar(i_out)/dvol(i_out)
mbar(i_out+1)            = m0(i_out+1)
r(i_out+1)               = r0(i_out+1)

!-----------------------------------------------------------------------
!              ||||| Radial zones i_out+1 to imax |||||
!-----------------------------------------------------------------------

DO i = i_out+1, imax
  dmbar(i)               = dm0(i)
  mbar(i+1)              = m0(i+1)
  vol(i+1)               = vol_0(i+1)
  dvol(i)                = dvol_0(i)
  rhobar(i)              = rhobar0(i)
  r(i+1)                 = x_ef(i+1)
END DO ! i = i_out+1, imax

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
!       RE-GRID TO CONSTANT DR/R FOR EACH SCALE HEIGHT IN DENSITY
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find the zones dividing density scale heights if the criteria for a
!   new regridding is satisfied. Otherwise, use the previously computed
!   grid.
!-----------------------------------------------------------------------

l_bounds                 = .false.
IF ( rhobar0(1) < rhobar_max ) THEN
  IF ( MOD( ncycle, int_pre_b  ) == 0 ) l_bounds = .true.
ELSE
  IF ( MOD( ncycle, int_post_b  ) == 0 ) l_bounds = .true.
END IF ! rhobar0 > rhobar_max

IF ( first ) THEN
  l_bounds               = .true.
  first                  = .false.
END IF

IF ( l_bounds ) THEN

  i_bounds               = 0
  i_10_hi                = EXPONENT( rhobar(i_in-1) )
  i_10_lo                = EXPONENT( rho_min        )
  n_scale                = i_10_hi - i_10_lo
  DO i = i_10_lo, i_10_hi
    ii                   = i_10_hi - i + 1
    i_bounds(ii)         = MAXLOC( rhobar(:), DIM = 1, MASK = EXPONENT(rhobar) == i )
  END DO !  i = i_10_lo, i_10_hi
  i_bounds(1)            = i_in
  i_bounds( i_10_hi-i_10_lo+2) = i_out

!-----------------------------------------------------------------------
!  If the lowest density partition contains 3 or fewer zones, combine
!   it with the next lowest partition
!-----------------------------------------------------------------------

  IF ( i_bounds(i_10_hi-i_10_lo+2) - i_bounds(i_10_hi-i_10_lo+1) <= 3 ) THEN
    i_bounds( i_10_hi-i_10_lo+1) = i_out
    n_scale              = n_scale - 1
  END IF

!-----------------------------------------------------------------------
!  If the highest density partition contains 3 or fewer zones, combine
!   it with the next highest partition
!-----------------------------------------------------------------------

  IF ( i_bounds(2) - i_bounds(1) <= 3 ) THEN
    DO i = 2, n_scale + 1
      i_bounds(i)        = i_bounds(i+1)
    END DO ! i = 3, n_scale + 1
    n_scale              = n_scale - 1
  END IF ! i_bounds(2) - i_bounds(1) <= 3
  
END IF ! MOD( ncycle, int_cycle ) == 0

DO i = 1, n_scale + 2
WRITE (nlog,3011) i, i_bounds(i)
 3011 FORMAT (' i=',i4,' i_bounds(i)=',i4)
END DO

!-----------------------------------------------------------------------
!  Smooth the grid in each density scale height
!-----------------------------------------------------------------------

DO i = 1, n_scale
  r_ratio                = ( r(i_bounds(i+1))/r(i_bounds(i)) )**( 1.d0/DBLE( i_bounds(i+1) - i_bounds(i) ) )
  DO ii = i_bounds(i) + 1, i_bounds(i+1)
    r(ii)                = r(ii-1) * r_ratio
  END DO ! ii = i_bounds(i) + 1, i_bounds(i+1)
END DO !  i = 1, n_scale + 1

!-----------------------------------------------------------------------
!  r(i+1) is set to x_ef(i+1), r(i) is calculated. Smooth the grid to
!   r(i+1) so that r(i) and r(i+1) remain adequately separated.
!-----------------------------------------------------------------------

i                        = n_scale + 1
r_ratio                  = ( r(i_bounds(i+1)+1)/r(i_bounds(i)) )**( 1.d0/DBLE( i_bounds(i+1) + 1 - i_bounds(i) ) )
DO ii = i_bounds(i) + 1, i_bounds(i+1) + 1
  r(ii)                  = r(ii-1) * r_ratio
END DO ! ii = i_bounds(i) + 1, i_bounds(i+1)

!-----------------------------------------------------------------------
!  Save rezoned radial grid
!-----------------------------------------------------------------------

r_new_grid(1:imax+1)     = r(1:imax+1)

!-----------------------------------------------------------------------
!  Limit change in r(i) to grid_frac * MIN( dr(i), dr(i-1) )
!-----------------------------------------------------------------------

l_repeat                 = .false.
DO i = 2, i_out
  IF ( DABS( r(i) - r0(i) ) > grid_frac * DMIN1( r0(i+1) - r0(i), r0(i) - r0(i-1) ) ) THEN
    r(i)                 = r0(i) + SIGN( grid_frac * DMIN1( r0(i+1) - r0(i), r0(i) - r0(i-1) ), r(i) - r0(i) )
    l_repeat             = .true.
  END IF ! r(i) - r0(i) > grid_frac * DMAX1( dr(i), dr(i-1) )
END DO ! 2, i_out

!-----------------------------------------------------------------------
!  Save rezoned radial grid after maximum change criteria imposed
!-----------------------------------------------------------------------

r0_prev(1:imax+1)        = r(1:imax+1)

!-----------------------------------------------------------------------
!  Compute new shell volumes and enclosed volumes
!-----------------------------------------------------------------------

DO i = 2, i_out - 1
  vol(i+1)               = frpith * r(i+1) * r(i+1) * r(i+1)
  dvol(i)                = vol(i+1) - vol(i)
END DO !  i = 2, i_out - 1
dvol(i_out)              = vol(i_out+1) - vol(i_out)

!-----------------------------------------------------------------------
!
!               \\\\\ NEW GRID HAS BEEN DETERMINED /////
!
!-----------------------------------------------------------------------

x_ef(imin:imax+1)        = r(imin:imax+1)
dx_cf(imin:imax)         = x_ef(imin+1:imax+1) - x_ef(imin:imax)
x_cf(imin:imax)          = half * ( x_ef(imin+1:imax+1) + x_ef(imin:imax) )

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
DEALLOCATE (vol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vol           '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol          '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (r0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r0            '; WRITE (nlog,2001) var_name; END IF
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

DEALLOCATE (rho_pad, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho_pad       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rhol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rhol          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rho6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rho6          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (drho, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'drho          '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dum, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dum           '; WRITE (nlog,2001) var_name; END IF

!-----------------------------------------------------------------------
!  Activate neutrino transport if high density zones have been adjusted
!-----------------------------------------------------------------------

  IF ( i_in /= i_in_d ) ncynu_trns = intnu_trns

RETURN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CONTAINS
REAL (KIND=double) FUNCTION rinterp( x1, x2, x3, y1, y2, y3, x )

REAL (KIND=double) :: x1
REAL (KIND=double) :: x2
REAL (KIND=double) :: x3
REAL (KIND=double) :: y1
REAL (KIND=double) :: y2
REAL (KIND=double) :: y3
REAL (KIND=double) :: x
REAL (KIND=double) :: a
REAL (KIND=double) :: b
REAL (KIND=double) :: c
REAL (KIND=double) :: denom

denom               = 1.d0/( ( x2 - x1 ) * ( x3 - x2 ) * ( x1 - x3 ) )

a                   = ( ( y2 - y1 ) * ( x3 - x1 ) - ( y3 - y1 ) * ( x2 - x1 ) ) * denom
b                   = ( ( y3 - y1 ) * ( x2 - x1 ) * ( x2 + x1 )           &
&                   -   ( y2 - y1 ) * ( x3 - x1 ) * ( x3 + x1 ) ) * denom
c                   = - ( y1 * x2 * x3 * ( x3 - x2 )                      &
&                   +     x1 * y2 * x3 * ( x1 - x3 )                      &
&                   +     x1 * x2 * y3 * ( x2 - x1 ) ) * denom
rinterp             = x * ( a * x + b ) + c
END FUNCTION rinterp

END SUBROUTINE regrid_final_grid

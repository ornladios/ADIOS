SUBROUTINE nuc_burn_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, ls, le, ldim, dgridp, tgridp, ygridp, rhoesp, idtyp, nprintp, &
& ncyclep, dtnph, time, rhop, tp, yep, xnp, a_nuc_repp, z_nuc_repp, &
& be_nuc_repp, nsep, dt, jdt )
!-----------------------------------------------------------------------
!
!    File:         nuc_burn_inout_fullnet
!    Module:       nuc_burn_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/16/04
!
!    Purpose:
!       To initialize the composition variables and the nuclear reaction
!        arrays
!
!    Subprograms called:
!  eqstt_x
!  network_setup
!  network_step
!  tgvndeye_comp_x
!
!    Input arguments:
!  imin          : lower x-array index
!  imax          : upper x-array index
!  nx            : x-array extent
!  ij_ray        : index denoting the j-index of a specific radial ray
!  ik_ray        : index denoting the k-index of a specific radial ray
!  ij_ray_dim    : number of y-zones on a processor before swapping
!  ik_ray_dim    : number of z-zones on a processor before swapping
!  ls            : lower composition index
!  le            : upper composition index
!  ldim          : composition array extent
!  dgridp        : number intervals in log rho per decade
!  tgridp        : number intervals in log t per decade
!  ygridp        : number intervals in ye between 0 and 1
!  rhoesp        : density separating eos gridding
!  idtyp         : gtid identifyer for radial zone j
!  nprintp       : unit number to print diagnostics
!  ncyclep       : cycle number
!  idtyp         : gtid identifyer for radial zone j
!  dtnph         : time step
!  rhop          : density (cm^{-3})
!  tp            : initial temprtature (K)
!  yep           : electron fraction
!  xnp           : initial mass fractions
!  a_nuc_repp    : mass number of the auxiliary heavy nucleus
!  z_nuc_repp    : charge number of the auxiliary heavy nucleus
!  be_nuc_repp   : binding energy of the auxiliary heavy nucleus
!  nesp          : nuclear statistical equilibrium flag
!
!    Output arguments:
!  tp            : updated temprtature (K)
!  xnp           : updated mass fractions
!  dt            : minimum nudlear burn time step restrictions (s)
!  jdt           : zones causing minimum time step restrictions dtp
!
!    Include files:
!  kind_module, numerical_module
!  incrmnt_module, nucbrn_module
!     
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one, epsilon

USE incrmnt_module, ONLY : dtmpmn
USE nucbrn_module, ONLY : dgrid, tgrid, ygrid, rhoes, idty, nprint, inuc, &
& dTmax, jdTmax, dynmax, jdynmax, dT_burn, dtime_burn, ncycle, dudt_nuc, &
& uburn, nse
USE eos_snc_x_module, ONLY : xn_e=>xn, a_nuc_rep, z_nuc_rep, be_nuc_rep

USE edit_module, ONLY : nprintt=>nprint


IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                  :: imin         ! lower x-array index
INTEGER, INTENT(in)                  :: imax         ! upper x-array index
INTEGER, INTENT(in)                  :: nx           ! x-array extent

INTEGER, INTENT(in)                  :: ij_ray       ! iindex denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                  :: ik_ray       ! iindex denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                  :: ij_ray_dim   ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                  :: ik_ray_dim   ! number of z-zones on a processor before swapping

INTEGER, INTENT(in)                  :: ls           ! inner logical composition index
INTEGER, INTENT(in)                  :: le           ! outer logical composition index
INTEGER, INTENT(in)                  :: ldim         ! physical dimension of x-zone

INTEGER, INTENT(in)                  :: ncyclep      ! cycle number
INTEGER, INTENT(in)                  :: nprintp      ! unit number to print diagnostics

INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: idtyp        ! index flag

REAL(KIND=double), INTENT(in), DIMENSION(3)              :: dgridp       ! eos density grid spacing
REAL(KIND=double), INTENT(in), DIMENSION(3)              :: tgridp       ! eos temperature grid spacing
REAL(KIND=double), INTENT(in), DIMENSION(3)              :: ygridp       ! eos electron fraction grid spacing
REAL(KIND=double), INTENT(in), DIMENSION(2)              :: rhoesp       ! density of grid changes
REAL(KIND=double), INTENT(in)                            :: dtnph        ! time step
REAL(KIND=double), INTENT(in)                            :: time         ! elapsed time
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: rhop         ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: yep          ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: a_nuc_repp   ! mass number of the auxiliary heavy nucleus
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: z_nuc_repp   ! charge number of the auxiliary heavy nucleus
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: be_nuc_repp  ! binding energy of the auxiliary heavy nucleus

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)           :: jdt          ! zone causing dtp

REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim) :: dt           ! minimum allowed time step

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)                :: nsep  ! nse flag

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: tp    ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ldim,ij_ray_dim,ik_ray_dim) :: xnp   ! mass fractions

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                               :: i               ! radial do index
INTEGER                               :: i_nc            ! shifted radial do index
INTEGER                               :: j               ! radial zone index
INTEGER                               :: l               ! composiiton index
INTEGER                               :: jr_min          ! minimum radial zone index (MGFLD indexing)
INTEGER                               :: jr_max          ! maximum radial zone index (MGFLD indexing)
INTEGER                               :: nonnse          ! index of first unshifted radial zone not in nse
INTEGER                               :: nuc_min         ! shifted index of first zone to compute nuclear burn
INTEGER                               :: nuc_max         ! shifted index of last zone to compute nuclear burn

INTEGER                               :: i_split         ! number of time splittings from tstart to tstop
INTEGER, PARAMETER                    :: i_split_max = 8 ! maximum number of time splittings to attempt from tstart to tstop
INTEGER                               :: n_subcycle      ! number of cycles to go from tstart to tstop
INTEGER                               :: i_subcycle      ! particular subcycle in going from tstart to tstop

REAL(KIND=double)                     :: tstart          ! time at beginning of time step
REAL(KIND=double)                     :: tstop           ! time at end of time step
REAL(KIND=double)                     :: esvd            ! dummy variable
REAL(KIND=double), DIMENSION(nx)      :: esvt            ! d(internal energy)/dT
REAL(KIND=double)                     :: esvy            ! dummy variable
REAL(KIND=double), DIMENSION(nx)      :: enm_initial     ! initial nuclear rest mass energy (ergs g^{-1})
REAL(KIND=double)                     :: enm_final       ! final nuclear rest mass energy (ergs g^{-1})
REAL(KIND=double)                     :: enb             ! nuclear binding energy (MeV/nucleon)

REAL(KIND=double)                     :: dtnph_sub       ! time etep of a subcycle
REAL(KIND=double)                     :: tstart_sub      ! time at the beginning of subcycle i_subcycle
REAL(KIND=double)                     :: tstop_sub       ! time at the end of subcycle i_subcycle

REAL(KIND=double), DIMENSION(nx)      :: rho             ! density (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)      :: t               ! temperature (K)
REAL(KIND=double), DIMENSION(nx)      :: t_i             ! initial temperature (K)
REAL(KIND=double), DIMENSION(nx)      :: ye              ! electron fraction
REAL(KIND=double), DIMENSION(nx)      :: e_i             ! internal energy (ergs g^&{-2})
REAL(KIND=double), DIMENSION(nx)      :: x_active        ! mass fraction of nuclei participating in reactions
REAL(KIND=double), DIMENSION(nx,ldim) :: xn              ! abundance mass fractions
REAL(KIND=double), DIMENSION(ldim)    :: xn_i            ! initial abundance mass fractions
REAL(KIND=double)                     :: x_sum           ! mass fraction sum
REAL(KIND=double)                     :: norm            ! normalization factor
REAL(KIND=double), PARAMETER          :: xn_tol = 0.1d0  ! maximum fractional change in xn allowed in a subcycle
REAL(KIND=double), PARAMETER          :: t_tol = 0.1d0   ! maximum fractional change in t allowed in a subcycle

 1001 FORMAT (' Error in closing nuclear_keys.d in subroutine nuc_burn_in')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                  \\\\\ SETUP FOR FULL_NET /////
!
!        Load variables received from radial_ray_module into
!         full_net modules.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize zone-center array size
!-----------------------------------------------------------------------

jr_min                      = 2
jr_max                      = imax - imin + 2

!-----------------------------------------------------------------------
!  Initialize time step parameters
!-----------------------------------------------------------------------

dTmax                       = zero
jdTmax                      = 0
dynmax                      = zero
jdynmax                     = 0
dtime_burn(1)               = 1.d+20
dtime_burn(2)               = 1.d+20

!-----------------------------------------------------------------------
!  Initialize temperature changes and energy generation rates
!-----------------------------------------------------------------------

dT_burn(jr_min:jr_max)      = zero
dudt_nuc(:,ij_ray,ik_ray)   = zero

!-----------------------------------------------------------------------
!  Return if inuc /= 1
!-----------------------------------------------------------------------

IF ( inuc /= 1 ) RETURN

!-----------------------------------------------------------------------
!  Cycle number
!-----------------------------------------------------------------------

ncycle                      = ncyclep

!-----------------------------------------------------------------------
!  Times
!-----------------------------------------------------------------------

tstart                      = time
tstop                       = time + dtnph

!-----------------------------------------------------------------------
!  Load nse flags in shifted array
!-----------------------------------------------------------------------

DO i = imin,imax
  nse(i-imin+2)             = nsep(i,ij_ray,ik_ray)
END DO

!-----------------------------------------------------------------------
!  Find nonnse, first ubshifted zone with nse = 0
!-----------------------------------------------------------------------

nonnse                      = 0
DO i = imin,imax
  IF ( nse(i+1) == 0 ) THEN
    nonnse                  = i
    EXIT
  END IF
END DO ! i

IF ( nonnse == 0 ) RETURN

!-----------------------------------------------------------------------
!        Indices for nuclear burning
!
!  Nuclear burn array indices are shifted from
!     nonnse : imax
!  to
!     1 : imax - nonnse + 1
!-----------------------------------------------------------------------

nuc_min                     = 1
nuc_max                     = imax - nonnse + 1

!-----------------------------------------------------------------------
!  Normalize and load abundances into MGFLD shifted array xn_e
!-----------------------------------------------------------------------

DO i_nc = nuc_min,nuc_max
  j                         = i_nc + nonnse
  x_sum                     = SUM(xnp(j-1,ls:le,ij_ray,ik_ray))
  norm                      = 1.d0/x_sum
  xnp(j-1,ls:le,ij_ray,ik_ray)  = norm * xnp(j-1,ls:le,ij_ray,ik_ray)
  xn_e(j,ls:le)             = xnp(j-1,ls:le,ij_ray,ik_ray)
END DO ! i_nc

!-----------------------------------------------------------------------
!  Load auxiliary heavy nucleus info into MGFLD shifted arrays
!-----------------------------------------------------------------------

DO i_nc = nuc_min,nuc_max
  j                         = i_nc + nonnse
  a_nuc_rep(j)              = a_nuc_repp (j-1,ij_ray,ik_ray)
  z_nuc_rep(j)              = z_nuc_repp (j-1,ij_ray,ik_ray)
  be_nuc_rep(j)             = be_nuc_repp(j-1,ij_ray,ik_ray)
END DO ! i_nc

!-----------------------------------------------------------------------
!  Load densities, temperatures and electron fractions, compute initial
!   internal energy and nuclear rest mass energy into nuc_burn shifted
!   arrays
!-----------------------------------------------------------------------

DO i = nonnse,imax
  i_nc                      = i-nonnse+1
  rho(i_nc)                 = rhop(i+imin-1,ij_ray,ik_ray)
  t(i_nc)                   = tp  (i+imin-1,ij_ray,ik_ray)
  t_i(i_nc)                 = t(i_nc)
  ye(i_nc)                  = yep (i+imin-1,ij_ray,ik_ray)
  j                         = i + 1
  CALL esrgn_comp_x( j, ij_ray, ik_ray, rho(i_nc), t(i_nc), ye(i_nc) )
  CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(i_nc), t(i_nc), ye(i_nc), &
&  e_i(i_nc), esvd, esvt(i_nc), esvy )
END DO ! i

!-----------------------------------------------------------------------
!  Load abundances into nuc_burn shifted arrays
!-----------------------------------------------------------------------

DO i = nonnse,imax
  xn(i-nonnse+1,ls:le)      = xnp(i+imin-1,ls:le,ij_ray,ik_ray)
  x_active(i-nonnse+1)      = one -  xn(i-nonnse+1,le) - xn(i-nonnse+1,le-1) &
&                           - xn(i-nonnse+1,le-2)
END DO ! i

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE FOR NUCLEAR BURN /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize temperature incrememts
!-----------------------------------------------------------------------

dT_burn(jr_min:jr_max)      = zero

!-----------------------------------------------------------------------
!
!           \\\\\ NUCLEAR BURN, SUBCYCLE IF NECESSARY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!            ||||| Begin loop over zones not in NSE |||||
!-----------------------------------------------------------------------

DO i = nonnse,imax
  i_nc                      = i-nonnse+1
  j                         = i + 1

!-----------------------------------------------------------------------
!  Store initial temperatures and mass fractions
!-----------------------------------------------------------------------

  t_i(i_nc)                 = tp(i+imin-1,ij_ray,ik_ray)
  xn_i(:)                   = xn(i_nc,:)

!-----------------------------------------------------------------------
!  Compute initial rest mass energy
!-----------------------------------------------------------------------

  CALL nuc_energy( j, enb, enm_initial(j) )

!-----------------------------------------------------------------------
!         ||||| Begin loop over time splitting attempts |||||
!-----------------------------------------------------------------------

  split : DO i_split = 0,i_split_max
    n_subcycle              = 2**i_split
    dtnph_sub               = dtnph/DBLE(n_subcycle)
    t(i_nc)                 = t_i(i_nc)
    xn(i_nc,:)              = xn_i(:)

!-----------------------------------------------------------------------
!                ||||| Begin loop over subcycles |||||
!-----------------------------------------------------------------------

    subcycle : DO i_subcycle = 1,n_subcycle
      tstart_sub            = tstart + dtnph_sub * DBLE( i_subcycle - 1 )
      tstop_sub             = tstart_sub + dtnph_sub

!-----------------------------------------------------------------------
!  Compute abundance change and energy generated by nuclear burn
!-----------------------------------------------------------------------

      CALL network_step( tstart_sub, tstop_sub, rho, t, ye, xn, i_nc, &
&      i_nc, nx, ldim, x_active, dtnph_sub )

      xn_e(j,ls:le)         = xn(i_nc,ls:le)

      CALL nuc_energy( j, enb, enm_final ) 
      t(i_nc)               = t_i(i_nc) - ( enm_final - enm_initial(j) )/esvt(i_nc)

!-----------------------------------------------------------------------
!  Test that abundance changes do not exceed allowed values
!-----------------------------------------------------------------------

      IF ( i_subcycle == 1 ) THEN
       IF ( DABS(t(i_nc) - t_i(i_nc))/( t(i_nc) ) > t_tol ) EXIT subcycle
        DO l = 1,ldim-3
          IF ( DABS(xn(i_nc,l) - xn_i(l))/( xn(i_nc,l) + 1.0d0 ) > xn_tol ) EXIT subcycle
        END DO ! l
      END IF ! i_subcycle == 1

      IF ( i_subcycle == n_subcycle ) EXIT split

!-----------------------------------------------------------------------
!                 ||||| End loop over subcycles |||||
!-----------------------------------------------------------------------

    END DO subcycle

!-----------------------------------------------------------------------
!          ||||| End loop over time splitting attempts |||||
!-----------------------------------------------------------------------

  END DO split

!-----------------------------------------------------------------------
!             ||||| End loop over zones not in NSE |||||
!-----------------------------------------------------------------------

END DO ! i

!-----------------------------------------------------------------------
!
!               \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Update temperature and compute nuclear energy released
!-----------------------------------------------------------------------

DO i_nc = nuc_min,nuc_max

  j                         = i_nc + nonnse

  DO l = ls,le
    xn_e(j,l)               = xn(i_nc,l)
  END DO

  CALL tgvndeye_comp_x( j, ij_ray, ik_ray, rho(i_nc), e_i(i_nc), ye(i_nc), &
&  t_i(i_nc), t(i_nc) )

  i                         = i_nc + nonnse - 1
  tp(i,ij_ray,ik_ray)       = t(i_nc)
  dT_burn(j)                = t(i_nc) - t_i(i_nc)
  dtmpmn(j,4,ij_ray,ik_ray) = dT_burn(j)

  CALL nuc_energy( j, enb, enm_final )
  uburn(j,ij_ray,ik_ray)    = uburn(j,ij_ray,ik_ray) - ( enm_final - enm_initial(j) )
  dudt_nuc(j,ij_ray,ik_ray) = - ( enm_final - enm_initial(j) )/( dtnph + epsilon )

END DO ! i_nc

!-----------------------------------------------------------------------
!  Updated abundances
!-----------------------------------------------------------------------

DO l = ls,le
  DO i_nc = nuc_min,nuc_max
    i                       = i_nc + nonnse - 1
    xnp(i,l,ij_ray,ik_ray)  = xn(i_nc,l)
  END DO
END DO

RETURN
END SUBROUTINE nuc_burn_inout

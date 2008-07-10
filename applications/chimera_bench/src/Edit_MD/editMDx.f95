SUBROUTINE editMDx( imin, imax, nx, jmin, jmax, ny, kmin, kmax, nz, nnu, &
& time, t_tb, y_in, z_in, r_nse, r_O1, r_xO, r_shock, r_shock_mn, r_shock_mx, &
& rsphere_mean, dsphere_mean, tsphere_mean, msphere_mean, esphere_mean, &
& tau_adv, tau_heat_nu, tau_heat_nuc, n_grow, r_gain, nprint )
!-----------------------------------------------------------------------
!
!    File:         editMDx
!    Module:       editMDx
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/22/05
!
!    Purpose:
!  To perform an edit of the entropies.
!
!    Subprograms call:
!  date_and_time_print : fetches and records the date and the time
!
!    Input arguments:
!  imin         : minimum x-array index for the edit
!  imax         : maximum x-array index for the edit
!  nx           : x-array extent
!  jmin         : minimum y-array index for the edit
!  jmax         : maximum y-array index for the edit
!  ny           : y_array extent
!  kmin         : minimum z-array index for the edit
!  kmax         : maximum z-array index for the edit
!  nz           : z dimension
!  nnu          : neutrino flavor array extent
!  time         : elapsed time (s)
!  t_tb         : time from core bounce (s)
!  y_in         : angular midpoint of zone
!  z_in         : z (azimuthal) midpoint of zone
!  r_nse        : radius of NSE-nonNSE boundary
!  r_O1         : radius of x(16O)=0.1 boundary
!  r_xO         : radius of x(16O)=0.5 boundary
!  r_shock      : radius of shock maximum
!  r_shock_mn   : minimum estimateed shock radius
!  r_shock_mx   : maximum estimateed shock radius
!  rsphere_mean : mean neutrinosphere radius (cm)
!  dsphere_mean : mean neutrinosphere density (g cm^{-3)}
!  tsphere_mean : mean neutrinosphere temperature (K)
!  msphere_mean : mean neutrinosphere enclosed mass (g)
!  esphere_mean : mean neutrinosphere rms neutrino energy (MeV)
!  tau_adv      : advection time scale (s)
!  tau_heat_nu  : neutrino heating time scale (s)
!  tau_heat_nuc : nuclear heating time scale (s) 
!  n_grow       : number of convective e-foldings in fluid from shock to gain radius
!  r_gain       : gain radius
!  iunitp       : unit number to print models
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                   :: imin           ! minimum x-array index for the edit
INTEGER, INTENT(in)                                   :: imax           ! maximum x-array index for the edit
INTEGER, INTENT(in)                                   :: nx             ! x-array extent

INTEGER, INTENT(in)                                   :: jmin           ! minimum y-array index for the edit
INTEGER, INTENT(in)                                   :: jmax           ! maximum y-array index for the edit
INTEGER, INTENT(in)                                   :: ny             ! y_array extent

INTEGER, INTENT(in)                                   :: kmin           ! minimum z-array index for the edit
INTEGER, INTENT(in)                                   :: kmax           ! maximum z-array index for the edit
INTEGER, INTENT(in)                                   :: nz             ! z_array extent

INTEGER, INTENT(in)                                   :: nnu            ! neutrino flavor array extent

INTEGER, INTENT(in)                                   :: nprint         ! unit number to print models

REAL(KIND=double), INTENT(in)                         :: time           ! elapsed time (s)
REAL(KIND=double), INTENT(in)                         :: t_tb           ! time from core bounce (s)
REAL(KIND=double), INTENT(in), DIMENSION(ny)          :: y_in           ! y (angular) midpoint of zone
REAL(KIND=double), INTENT(in), DIMENSION(nz)          :: z_in           ! z (azimuthal) midpoint of zone

REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: r_nse          ! radius of NSE-nonNSE boundary
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: r_O1           ! radius of x(16O)=0.1 boundary
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: r_xO           ! radius of x(16O)=0.5 boundary
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: r_shock        ! radius of shock maximum
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: r_shock_mn     ! minimum estimateed shock radius
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: r_shock_mx     ! maximum estimateed shock radius
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)   :: rsphere_mean   ! mean neutrinosphere radius
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)   :: dsphere_mean   ! mean neutrinosphere density
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)   :: tsphere_mean   ! mean neutrinosphere temperature
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)   :: msphere_mean   ! mean neutrinosphere enclosed mass
REAL(KIND=double), INTENT(in), DIMENSION(nnu,ny,nz)   :: esphere_mean   ! mean neutrinosphere energy
REAL(KIND=double), INTENT(in), DIMENSION(nnu+1,ny,nz) :: r_gain         ! gain radius
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: tau_adv        ! advection time scale (s)
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: tau_heat_nu    ! neutrino heating time scale (s)
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: tau_heat_nuc   ! nuclear heating time scale (s) 
REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)       :: n_grow         ! number of convective e-foldings in fluid from shock to gain radius

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=11)                                    :: name             ! place holder for name of variable

INTEGER, PARAMETER                                    :: i = 0            ! integer placement holder

REAL(KIND=double)                                     :: r_nse_min        ! minimum of r_nse
REAL(KIND=double)                                     :: r_nse_max        ! maximum of r_nse
REAL(KIND=double)                                     :: r_O1_min         ! minimum of r_O1
REAL(KIND=double)                                     :: r_O1_max         ! maximum of r_O1
REAL(KIND=double)                                     :: r_xO_min         ! minimum of r_O5
REAL(KIND=double)                                     :: r_xO_max         ! maximum of r_O5
REAL(KIND=double)                                     :: r_shk_min        ! minimum of r_shock
REAL(KIND=double)                                     :: r_shk_max        ! maximum of r_shock
REAL(KIND=double)                                     :: r_shk_mn_min     ! minimum of r_shock_mn
REAL(KIND=double)                                     :: r_shk_mn_max     ! maximum of r_shock_mn
REAL(KIND=double)                                     :: r_shk_mx_min     ! minimum of r_shock_mx
REAL(KIND=double)                                     :: r_shk_mx_max     ! maximum of r_shock_mx
REAL(KIND=double)                                     :: rsphere_min      ! minimum of rsphere_mean
REAL(KIND=double)                                     :: rsphere_max      ! maximum of rsphere_mean
REAL(KIND=double)                                     :: dsphere_min      ! minimum of dsphere_mean
REAL(KIND=double)                                     :: dsphere_max      ! maximum of dsphere_mean
REAL(KIND=double)                                     :: tsphere_min      ! minimum of tsphere_mean
REAL(KIND=double)                                     :: tsphere_max      ! maximum of tsphere_mean
REAL(KIND=double)                                     :: msphere_min      ! minimum of msphere_mean
REAL(KIND=double)                                     :: msphere_max      ! maximum of msphere_mean
REAL(KIND=double)                                     :: esphere_min      ! minimum of esphere_mean
REAL(KIND=double)                                     :: esphere_max      ! maximum of esphere_mean
REAL(KIND=double)                                     :: r_gain_min       ! minimum of r_gain
REAL(KIND=double)                                     :: r_gain_max       ! maximum of r_gain
REAL(KIND=double)                                     :: tau_adv_min      ! minimum of tau_adv
REAL(KIND=double)                                     :: tau_adv_max      ! maximum of tau_adv
REAL(KIND=double)                                     :: tau_heat_nu_min  ! minimum of tau_heat_nu
REAL(KIND=double)                                     :: tau_heat_nu_max  ! maximum of tau_heat_nu
REAL(KIND=double)                                     :: tau_heat_nuc_min ! minimum of tau_heat_nuc
REAL(KIND=double)                                     :: tau_heat_nuc_max ! maximum of tau_heat_nuc
REAL(KIND=double)                                     :: n_grow_min       ! minimum of tau_heat_nuc
REAL(KIND=double)                                     :: n_grow_max       ! maximum of tau_heat_nuc

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

WRITE (nprint) imin, imax, nx, jmin, jmax, ny
WRITE (nprint) time, t_tb
WRITE (nprint) y_in

name              = '  r_nse    '
r_nse_min         = MINVAL( r_nse(:,:) )
r_nse_max         = MAXVAL( r_nse(:,:) )
WRITE (nprint) i,name,r_nse(:,:),r_nse_min,r_nse_max

name              = '  r_O1     '
r_O1_min          = MINVAL( r_O1(:,:) )
r_O1_max          = MAXVAL( r_O1(:,:) )
WRITE (nprint) i,name,r_O1(:,:),r_O1_min,r_O1_max

name              = '  r_O5     '
r_xO_min          = MINVAL( r_xO(:,:) )
r_xO_max          = MAXVAL( r_xO(:,:) )
WRITE (nprint) i,name,r_xO(:,:),r_xO_min,r_xO_max

name              = ' r_shock   '
r_shk_min         = MINVAL( r_shock(:,:) )
r_shk_max         = MAXVAL( r_shock(:,:) )
WRITE (nprint) i,name,r_shock(:,:),r_shk_min,r_shk_max

name              = ' r_shock_mn'
r_shk_mn_min      = MINVAL( r_shock_mn(:,:) )
r_shk_mn_max      = MAXVAL( r_shock_mn(:,:) )
WRITE (nprint) i,name,r_shock_mn(:,:),r_shk_mn_min,r_shk_mn_max

name              = ' r_shock_mx'
r_shk_mx_min      = MINVAL( r_shock_mx(:,:) )
r_shk_mx_max      = MAXVAL( r_shock_mx(:,:) )
WRITE (nprint) i,name,r_shock_mx(:,:),r_shk_mx_min,r_shk_mx_max

name              = ' r_sphere_1'
rsphere_min       = MINVAL( rsphere_mean(1,:,:) )
rsphere_max       = MAXVAL( rsphere_mean(1,:,:) )
WRITE (nprint) i,name,rsphere_mean(1,:,:),rsphere_min,rsphere_max

name              = ' r_sphere_2'
rsphere_min       = MINVAL( rsphere_mean(2,:,:) )
rsphere_max       = MAXVAL( rsphere_mean(2,:,:) )
WRITE (nprint) i,name,rsphere_mean(2,:,:),rsphere_min,rsphere_max

name              = ' r_sphere_3'
rsphere_min       = MINVAL( rsphere_mean(3,:,:) )
rsphere_max       = MAXVAL( rsphere_mean(3,:,:) )
WRITE (nprint) i,name,rsphere_mean(3,:,:),rsphere_min,rsphere_max

name              = ' r_sphere_4'
rsphere_min       = MINVAL( rsphere_mean(4,:,:) )
rsphere_max       = MAXVAL( rsphere_mean(4,:,:) )
WRITE (nprint) i,name,rsphere_mean(4,:,:),rsphere_min,rsphere_max

name              = ' d_sphere_1'
dsphere_min       = MINVAL( dsphere_mean(1,:,:) )
dsphere_max       = MAXVAL( dsphere_mean(1,:,:) )
WRITE (nprint) i,name,dsphere_mean(1,:,:),dsphere_min,dsphere_max

name              = ' d_sphere_2'
dsphere_min       = MINVAL( dsphere_mean(2,:,:) )
dsphere_max       = MAXVAL( dsphere_mean(2,:,:) )
WRITE (nprint) i,name,dsphere_mean(2,:,:),dsphere_min,dsphere_max

name              = ' d_sphere_3'
dsphere_min       = MINVAL( dsphere_mean(3,:,:) )
dsphere_max       = MAXVAL( dsphere_mean(3,:,:) )
WRITE (nprint) i,name,dsphere_mean(3,:,:),dsphere_min,dsphere_max

name              = ' d_sphere_4'
dsphere_min       = MINVAL( dsphere_mean(4,:,:) )
dsphere_max       = MAXVAL( dsphere_mean(4,:,:) )
WRITE (nprint) i,name,dsphere_mean(4,:,:),dsphere_min,dsphere_max

name              = ' t_sphere_1'
tsphere_min       = MINVAL( tsphere_mean(1,:,:) )
tsphere_max       = MAXVAL( tsphere_mean(1,:,:) )
WRITE (nprint) i,name,tsphere_mean(1,:,:),tsphere_min,tsphere_max

name              = ' t_sphere_2'
tsphere_min       = MINVAL( tsphere_mean(2,:,:) )
tsphere_max       = MAXVAL( tsphere_mean(2,:,:) )
WRITE (nprint) i,name,tsphere_mean(2,:,:),tsphere_min,tsphere_max

name              = ' t_sphere_3'
tsphere_min       = MINVAL( tsphere_mean(3,:,:) )
tsphere_max       = MAXVAL( tsphere_mean(3,:,:) )
WRITE (nprint) i,name,tsphere_mean(3,:,:),tsphere_min,tsphere_max

name              = ' t_sphere_4'
tsphere_min       = MINVAL( tsphere_mean(4,:,:) )
tsphere_max       = MAXVAL( tsphere_mean(4,:,:) )
WRITE (nprint) i,name,tsphere_mean(4,:,:),tsphere_min,tsphere_max

name              = ' m_sphere_1'
msphere_min       = MINVAL( msphere_mean(1,:,:) )
msphere_max       = MAXVAL( msphere_mean(1,:,:) )
WRITE (nprint) i,name,msphere_mean(1,:,:),msphere_min,msphere_max

name              = ' m_sphere_2'
msphere_min       = MINVAL( msphere_mean(2,:,:) )
msphere_max       = MAXVAL( msphere_mean(2,:,:) )
WRITE (nprint) i,name,msphere_mean(2,:,:),msphere_min,msphere_max

name              = ' m_sphere_3'
msphere_min       = MINVAL( msphere_mean(3,:,:) )
msphere_max       = MAXVAL( msphere_mean(3,:,:) )
WRITE (nprint) i,name,msphere_mean(3,:,:),msphere_min,msphere_max

name              = ' m_sphere_4'
msphere_min       = MINVAL( msphere_mean(4,:,:) )
msphere_max       = MAXVAL( msphere_mean(4,:,:) )
WRITE (nprint) i,name,msphere_mean(4,:,:),msphere_min,msphere_max

name              = ' e_sphere_1'
esphere_min       = MINVAL( esphere_mean(1,:,:) )
esphere_max       = MAXVAL( esphere_mean(1,:,:) )
WRITE (nprint) i,name,esphere_mean(1,:,:),esphere_min,esphere_max

name              = ' e_sphere_2'
esphere_min       = MINVAL( esphere_mean(2,:,:) )
esphere_max       = MAXVAL( esphere_mean(2,:,:) )
WRITE (nprint) i,name,esphere_mean(2,:,:),esphere_min,esphere_max

name              = ' e_sphere_3'
esphere_min       = MINVAL( esphere_mean(3,:,:) )
esphere_max       = MAXVAL( esphere_mean(3,:,:) )
WRITE (nprint) i,name,esphere_mean(3,:,:),esphere_min,esphere_max

name              = ' e_sphere_4'
esphere_min       = MINVAL( esphere_mean(4,:,:) )
esphere_max       = MAXVAL( esphere_mean(4,:,:) )
WRITE (nprint) i,name,esphere_mean(4,:,:),esphere_min,esphere_max

name              = '  r_gain_1 '
r_gain_min        = MINVAL( r_gain(1,:,:) )
r_gain_max        = MAXVAL( r_gain(1,:,:) )
WRITE (nprint) i,name,r_gain(1,:,:),r_gain_min,r_gain_max

name              = '  r_gain_2 '
r_gain_min        = MINVAL( r_gain(2,:,:) )
r_gain_max        = MAXVAL( r_gain(2,:,:) )
WRITE (nprint) i,name,r_gain(2,:,:),r_gain_min,r_gain_max

name              = '  r_gain   '
r_gain_min        = MINVAL( r_gain(5,:,:) )
r_gain_max        = MAXVAL( r_gain(5,:,:) )
WRITE (nprint) i,name,r_gain(5,:,:),r_gain_min,r_gain_max

name              = '  tau_adv  '
tau_adv_min       = MINVAL( tau_adv(:,:) )
tau_adv_max       = MAXVAL( tau_adv(:,:) )
WRITE (nprint) i,name,tau_adv(:,:),tau_adv_min,tau_adv_max

name              = ' tau_q_nu  '
tau_heat_nu_min   = MINVAL( tau_heat_nu(:,:) )
tau_heat_nu_max   = MAXVAL( tau_heat_nu(:,:) )
WRITE (nprint) i,name,tau_heat_nu(:,:),tau_heat_nu_min,tau_heat_nu_max

name              = ' tau_q_nuc '
tau_heat_nuc_min  = MINVAL( tau_heat_nuc(:,:) )
tau_heat_nuc_max  = MAXVAL( tau_heat_nuc(:,:) )
WRITE (nprint) i,name,tau_heat_nuc(:,:),tau_heat_nuc_min,tau_heat_nuc_max

name              = '  cv_grow  '
n_grow_min        = MINVAL( n_grow(:,:) )
n_grow_max        = MAXVAL( n_grow(:,:) )
WRITE (nprint) i,name,n_grow(:,:),n_grow_min,n_grow_max

RETURN
END  SUBROUTINE editMDx

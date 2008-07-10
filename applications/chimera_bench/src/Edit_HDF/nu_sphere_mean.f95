SUBROUTINE nu_sphere_mean( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, e_rms_trns, j_sphere, r_sphere, d_sphere, t_sphere, m_sphere, &
& nx, nez, nnu, rsphere_mean, dsphere_mean, tsphere_mean, msphere_mean, &
& esphere_mean, jsphere_mean )
!-----------------------------------------------------------------------
!
!    File:         nu_sphere_mean
!    Module:       nu_sphere_mean
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/22/00
!
!    Purpose:
!      To calculate the radius, density, temperature, and enclosed
!       mass of the neutrino-spheres evaluated at the rms neutrino
!       energy.
!
!      Subroutine nu_sphere and e_rms must be called first.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min                : inner radial zone number
!  jr_max                : outer radial zone number
!  ij_ray                : j-index of a radial ray
!  ik_ray                : k-index of a radial ray
!  ij_ray_dim            : number of y-zones on a processor before swapping with y
!  ik_ray_dim            : number of z-zones on a processor before swapping with z
!  e_rms_trns            : rms transport neutrino energy
!  j_sphere              : radial indices of the (k,n) neutrinospheres
!  r_sphere              : radii of the (k,n) neutrinospheres
!  d_sphere              : densities of the (k,n) neutrinospheres
!  t_sphere              : temperatures of the (k,n) neutrinospheres
!  m_sphere              : enclosed masses of the (k,n) neutrinospheres
!  nx                    : x_array extent
!  nez                   : neutrino energy array extent
!  nnu                   : neutrino flavor array extent
!
!    Output arguments:
!  jsphere_mean(n,ij_ray,ik_ray) : shifted radial zone <= mean n-neutrinosphere
!  rsphere_mean(n,ij_ray,ik_ray) : radius of the mean n-neutrinosphere
!  dsphere_mean(n,ij_ray,ik_ray) : density of the mean n-neutrinosphere
!  tsphere_mean(n,ij_ray,ik_ray) : temperature of the mean n-neutrinosphere
!  msphere_mean(n,ij_ray,ik_ray) : enclosed rest mass of the mean n-neutrinosphere
!
!    Output arguments (common):
!
!    Include files:
!  kind_module, array_module, numerical_module
!  mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half, one, epsilon

USE nu_dist_module, ONLY : unu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: nx            ! radial array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

INTEGER, INTENT(in), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim)           :: j_sphere   ! zone index of the neutrinospheres

REAL(KIND=double), INTENT(in), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: r_sphere   ! neutrinosphere radius
REAL(KIND=double), INTENT(in), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: d_sphere   ! neutrinosphere density
REAL(KIND=double), INTENT(in), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: t_sphere   ! neutrinosphere temperature
REAL(KIND=double), INTENT(in), DIMENSION(nez,nnu,ij_ray_dim,ik_ray_dim) :: m_sphere   ! neutrinosphere enclosed mass
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)  :: e_rms_trns ! SUM psi1 * w5dw/SUM w3ww

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)           :: jsphere_mean  ! mean neutrinosphere radius

REAL(KIND=double), INTENT(out), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: rsphere_mean  ! mean neutrinosphere radius
REAL(KIND=double), INTENT(out), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: dsphere_mean  ! mean neutrinosphere density
REAL(KIND=double), INTENT(out), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: tsphere_mean  ! mean neutrinosphere temperature
REAL(KIND=double), INTENT(out), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: msphere_mean  ! mean neutrinosphere enclosed mass
REAL(KIND=double), INTENT(out), DIMENSION(nnu,ij_ray_dim,ik_ray_dim) :: esphere_mean  ! mean neutrinosphere energy

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index

INTEGER                          :: k_rms         ! smallest neutrino energy index where e_rms < e_sphere
INTEGER                          :: k_rms_1       ! k_rms - 1
INTEGER                          :: j_sph         ! j_sphere(k_rms,n,ij_ray,ik_ray)
INTEGER                          :: j_sph_1       ! j_sphere(k_rms_1,n,ij_ray,ik_ray)

REAL(KIND=double)                :: e_sphere      ! neutrino energy at j_sphere
REAL(KIND=double)                :: e_sphere_1    ! neutrino energy at j_sphere
REAL(KIND=double)                :: e_rms         ! e_rms_trns(j_sph,n)
REAL(KIND=double)                :: e_rms_1       ! e_rms_trns(j_sph_1,n)

REAL(KIND=double)                :: e_diff        ! e_sphere - e_rms
REAL(KIND=double)                :: e_diff_1      ! e_sphere_1 - e_rms_1
REAL(KIND=double)                :: r_intrp       ! interpolation coefficient

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------

jsphere_mean(:,ij_ray,ik_ray) = 1
rsphere_mean(:,ij_ray,ik_ray) = zero
dsphere_mean(:,ij_ray,ik_ray) = zero
tsphere_mean(:,ij_ray,ik_ray) = zero
msphere_mean(:,ij_ray,ik_ray) = zero
esphere_mean(:,ij_ray,ik_ray) = zero

!-----------------------------------------------------------------------
!        Return if no neutrinos.
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!
!           \\\\\ MEAN N-NEUTRINO EMISSION SUEFACES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Find energy index k_rms such that
!
!    unu(j_sph,k_rms-1) <  e_rms_trns(j_sph,n,ij_ray,ik_ray) <  unu(j_sph,k_rms)
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k = 1,nnugp(n)
    k_rms                   = k
    j_sph                   = j_sphere(k,n,ij_ray,ik_ray)
    e_sphere                = unu(j_sph,k)
    e_rms                   = e_rms_trns(j_sph,n,ij_ray,ik_ray)
    IF ( e_rms < e_sphere ) EXIT
  END DO ! k = 1,nnugp(n)

  k_rms_1                   = MAX( k_rms - 1, 1 )  
    
  IF ( k_rms == k_rms_1 ) THEN
    jsphere_mean(n,ij_ray,ik_ray)   = j_sph_1
    rsphere_mean(n,ij_ray,ik_ray)   = r_sphere(k_rms,n,ij_ray,ik_ray)
    dsphere_mean(n,ij_ray,ik_ray)   = d_sphere(k_rms,n,ij_ray,ik_ray)
    tsphere_mean(n,ij_ray,ik_ray)   = t_sphere(k_rms,n,ij_ray,ik_ray)
    msphere_mean(n,ij_ray,ik_ray)   = m_sphere(k_rms,n,ij_ray,ik_ray)
    esphere_mean(n,ij_ray,ik_ray)   = e_sphere
  ELSE ! k_rms /= k_rms_1
    j_sph_1                 = j_sphere(k_rms_1,n,ij_ray,ik_ray)
    e_sphere_1              = unu(j_sph,k_rms_1)
    e_rms_1                 = e_rms_trns(j_sph_1,n,ij_ray,ik_ray)
    e_diff                  = e_sphere   - e_rms
    e_diff_1                = e_sphere_1 - e_rms_1
    r_intrp                 = e_diff/( e_diff - e_diff_1 )

    jsphere_mean(n,ij_ray,ik_ray)   = j_sph_1
    rsphere_mean(n,ij_ray,ik_ray)   = r_sphere(k_rms,n,ij_ray,ik_ray)  + ( r_sphere(k_rms_1,n,ij_ray,ik_ray)  &
&                                   - r_sphere(k_rms,n,ij_ray,ik_ray) ) * r_intrp
    dsphere_mean(n,ij_ray,ik_ray)   = d_sphere(k_rms,n,ij_ray,ik_ray)  + ( d_sphere(k_rms_1,n,ij_ray,ik_ray)  &
&                                   - d_sphere(k_rms,n,ij_ray,ik_ray) ) * r_intrp
    tsphere_mean(n,ij_ray,ik_ray)   = t_sphere(k_rms,n,ij_ray,ik_ray)  + ( t_sphere(k_rms_1,n,ij_ray,ik_ray)  &
&                                   - t_sphere(k_rms,n,ij_ray,ik_ray) ) * r_intrp
    msphere_mean(n,ij_ray,ik_ray)   = m_sphere(k_rms,n,ij_ray,ik_ray)  + ( m_sphere(k_rms_1,n,ij_ray,ik_ray)  &
&                                   - m_sphere(k_rms,n,ij_ray,ik_ray) ) * r_intrp
    esphere_mean(n,ij_ray,ik_ray)   = e_sphere                         + ( e_sphere_1                         &
&                                   - e_sphere )                        * r_intrp
  END IF ! k_rms == k_rms_1

END DO ! n = 1,nnu

RETURN
END SUBROUTINE nu_sphere_mean

SUBROUTINE Brunt_Vaisala_yl( j, ij_ray, ik_ray, nx, r, rho, yl, g_acc_e, &
& wBV, twBV, wBV_s, twBV_s, wBV_yl, twBV_yl )
!-----------------------------------------------------------------------
!
!    File:         Brunt_Vaisala_yl
!    Module:       Brunt_Vaisala_yl
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/01/08
!
!    Purpose:
!      To compute the the Brunt-Vaisala angular frequency between
!       adjacent radial zones j and j+2 as a function of gradiants of
!       s and yl (i.e., ds/dr and dyl/dr) for the neutrino opaque
!       case. In this case the matter and neutrino fluids are treated
!       as one fluid, and the equation of state and composition is taken
!       to be a function of p, s, and yl.
!
!    Input arguments:
!  j         : shifted radial zone index
!  ij_ray    : j-index of a radial ray
!  ik_ray    : k-index of a radial ray
!  nx        : x-array extent
!  r         : zone-edged radii [cm]
!  rho       : density [g cm^{-3}]
!  yl        : lepton fraction
!  g_acc_e   : zone-edged x-component of gravitational acceleration [cm s^{-2} g^{-1}]
!
!    Output arguments:
!  wBV       : Brunt-Vaisala frequency (positive if unstable)
!  twBV      : Brunt-Vaisala growth time (positive if unstable)
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  eos_drv_module, eos_snc_x_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half, epsilon

USE eos_drv_module, ONLY : aesv, aesvd, aesvt, aesvy
USE eos_snc_x_module, ONLY : nse

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)            :: j            ! shifted radial zone index
INTEGER, INTENT(in)            :: ij_ray       ! j-index of a radial ray
INTEGER, INTENT(in)            :: ik_ray       ! k-index of a radial ray
INTEGER, INTENT(in)            :: nx           ! x-array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: r        ! zone-edged radii [cm]
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho      ! density [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: yl       ! lepton fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: g_acc_e  ! zone-edged x-component of gravitational acceleration (cm s^{-2} g^{-1})

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out) :: wBV          ! Brunt-Vaisala frequency (positive if unstable) [s^{-1}]
REAL(KIND=double), INTENT(out) :: twBV         ! Brunt-Vaisala growth time (positive if unstable) [s]
REAL(KIND=double), INTENT(out) :: wBV_s        ! Brunt-Vaisala frequency due to entropy (positive if unstable) [s^{-1}]
REAL(KIND=double), INTENT(out) :: twBV_s       ! Brunt-Vaisala growth time due to entropy (positive if unstable) [s]
REAL(KIND=double), INTENT(out) :: wBV_yl       ! Brunt-Vaisala frequency due to lepton fraction (positive if unstable) [s^{-1}]
REAL(KIND=double), INTENT(out) :: twBV_yl      ! Brunt-Vaisala growth time due to lepton fraction (positive if unstable) [s]

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)              :: dpdd         ! d(presure)/(density)
REAL(KIND=double)              :: dpdt         ! d(presure)/(temperature)
REAL(KIND=double)              :: dpdy         ! d(presure)/(electron fraction)
REAL(KIND=double)              :: dsdd         ! d(entropy)/(density)
REAL(KIND=double)              :: dsdt         ! d(entropy)/(temperature)
REAL(KIND=double)              :: dsdy         ! d(entropy)/(electron fraction)

REAL(KIND=double)              :: denom        ! denominator in Cramer's Law
REAL(KIND=double)              :: drodyl_ps    ! d(density)/d(electron fraction)
REAL(KIND=double)              :: drods_pyl    ! d(density)/d(entropy)
REAL(KIND=double)              :: dsdr         ! d(entropy)/d(radius)
REAL(KIND=double)              :: dyldr        ! d(electron fraction)/d(radius)

REAL(KIND=double)              :: sj           ! enropy at radial zone j
REAL(KIND=double)              :: sjp1         ! enropy at radial zone j+1

REAL(KIND=double)              :: w2BV         ! square of the Brunt-Vaisala frequency
REAL(KIND=double)              :: wBV2_s       ! square of the Brunt-Vaisala frequency
REAL(KIND=double)              :: wBV2_yl      ! square of the Brunt-Vaisala frequency

REAL(KIND=double)              :: rho_e        ! zone-edged density [g cm^{-3}]
REAL(KIND=double)              :: r_c_j        ! zone-centered radii at j [cm]
REAL(KIND=double)              :: r_c_jp1      ! zone-centered radii at j+1 [cm]

!-----------------------------------------------------------------------
!  Zero quantities and return if nse = 0
!-----------------------------------------------------------------------

wBV                = zero
twBV               = zero
wBV_s              = zero
twBV_s             = zero
wBV_yl             = zero
twBV_yl            = zero

IF ( nse(j,ij_ray,ik_ray) == 0 ) RETURN

!-----------------------------------------------------------------------
!  Zone-centered radii
!-----------------------------------------------------------------------

r_c_j              = half * ( r(j  ) + r(j-1) )
r_c_jp1            = half * ( r(j+1) + r(j  ) )

!-----------------------------------------------------------------------
!  Zone-edged density
!-----------------------------------------------------------------------

rho_e              = half * ( rho(j) + rho(j+1) )

!-----------------------------------------------------------------------
!  Zone edged entropies and derivatives
!-----------------------------------------------------------------------

dpdd               = half * ( aesvd(j,1,ij_ray,ik_ray) + aesvd(j+1,1,ij_ray,ik_ray) )
dpdt               = half * ( aesvt(j,1,ij_ray,ik_ray) + aesvt(j+1,1,ij_ray,ik_ray) )
dpdy               = half * ( aesvy(j,1,ij_ray,ik_ray) + aesvy(j+1,1,ij_ray,ik_ray) )
sj                 = aesv (j  ,3,ij_ray,ik_ray)
sjp1               = aesv (j+1,3,ij_ray,ik_ray)
dsdd               = half * ( aesvd(j,3,ij_ray,ik_ray) + aesvd(j+1,3,ij_ray,ik_ray) )
dsdt               = half * ( aesvt(j,3,ij_ray,ik_ray) + aesvt(j+1,3,ij_ray,ik_ray) )
dsdy               = half * ( aesvy(j,3,ij_ray,ik_ray) + aesvy(j+1,3,ij_ray,ik_ray) )

!-----------------------------------------------------------------------
!  Compute d(rho)/d(s) and d(rho)/d(ye)
!-----------------------------------------------------------------------
!
!                       s  drho |
!        drods_pyl  :  --- ---- |
!                      rho  ds  |ye,p
!
!                      ye  drho |
!        drodyl_ps  :  --- ---- |
!                      rho dye  |s,p
!
!-----------------------------------------------------------------------

denom              = dsdd * dpdt - dsdt * dpdd
IF ( DABS(denom) < epsilon ) RETURN 
drodyl_ps          = ( dpdy * dsdt - dpdt * dsdy )/denom
drods_pyl          = dpdt/denom

!-----------------------------------------------------------------------
!  Compute the radial gradiants dsdr and dyldr
!-----------------------------------------------------------------------

dsdr               = ( sjp1    - sj    )/( r_c_jp1 - r_c_j )
dyldr              = ( yl(j+1) - yl(j) )/( r_c_jp1 - r_c_j )

!-----------------------------------------------------------------------
!  Compute the Brunt-Vaisala frequency and time scale
!-----------------------------------------------------------------------

w2BV               = - ( g_acc_e(j)/rho_e )                             &
&                  * ( drods_pyl * dsdr + drodyl_ps * dyldr )
wBV                = SIGN( DSQRT( DABS(w2BV) + epsilon ), w2BV )
twBV               = 1.d+0/( wBV + epsilon )
wBV2_s             = - ( g_acc_e(j)/rho_e ) * ( drods_pyl * dsdr )
wBV_s              = SIGN( DSQRT( DABS(wBV2_s) + epsilon ), wBV2_s )
twBV_s             = 1.d+0/( wBV_s + epsilon )
wBV2_yl            = - ( g_acc_e(j)/rho_e ) * ( drodyl_ps * dyldr )
wBV_yl             = SIGN( DSQRT( DABS(wBV2_yl) + epsilon ), wBV2_yl )
twBV_yl            = 1.d+0/( wBV_yl + epsilon )

RETURN
END SUBROUTINE Brunt_Vaisala_yl


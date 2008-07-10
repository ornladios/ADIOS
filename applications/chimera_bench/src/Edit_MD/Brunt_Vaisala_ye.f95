SUBROUTINE Brunt_Vaisala_ye( i, ij_ray, ik_ray, nx, ij_ray_dim, ik_ray_dim, &
& x_c, rho_c, ye_c, grav_x_c, wBV, twBV )
!-----------------------------------------------------------------------
!
!    File:         Brunt-Vaisala_ye
!    Module:       Brunt-Vaisala_ye
!    Type:         program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/11/96
!
!    Purpose:
!      To compute the the Brunt-Vaisala angular frequency as a
!       function of gradiants of s and ye (i.e., ds/dr and dye/dr)
!       for the neutrino transpoarent case. In this case the matter
!       and neutrino fluids are treated as one fluid, and the equation
!       of state and composition is taken to be a function of p,
!       s, and ye.
!
!    Input arguments:
!  i          : unshifted radial zone index
!  i_ray      : index denoting a specific radial ray
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!  nx         : x-array extent
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  x_c        : radial midpoint of zone [cm]
!  rho_c      : density of zone [g cm^{-3}]
!  ye_c       : electron fraction
!  grav_x_c   : zone-centered x-component of gravitational acceleration [cm s^{-2} g^{-1}]
!
!    Output arguments:
!  wBV        : Brunt-Vaisala frequency (positive if unstable)
!  twBV       : Brunt-Vaisala growth time (positive if unstable)
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module
!  eos_snc_x_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE eos_snc_x_module, ONLY : aesv, aesvd, aesvt, aesvy

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)            :: i            ! unshifted radial zone index
INTEGER, INTENT(in)            :: ij_ray       ! j-index of a radial ray
INTEGER, INTENT(in)            :: ik_ray       ! k-index of a radial ray
INTEGER, INTENT(in)            :: nx           ! x-array extent
INTEGER, INTENT(in)            :: ij_ray_dim   ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)            :: ik_ray_dim   ! number of z-zones on a processor before swapping with z

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: x_c             ! radial midpoint of zone [cm]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: rho_c    ! density of zone [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: ye_c     ! entropy of zone
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: grav_x_c ! zone-centered x-component of gravitational acceleration (cm s^{-2} g^{-1})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                        :: j            ! shifted radial zone index

REAL(KIND=double)              :: dpdd         ! d(presure)/(density)
REAL(KIND=double)              :: dpdt         ! d(presure)/(temperature)
REAL(KIND=double)              :: dpdy         ! d(presure)/(electron fraction)
REAL(KIND=double)              :: dsdd         ! d(entropy)/(density)
REAL(KIND=double)              :: dsdt         ! d(entropy)/(temperature)
REAL(KIND=double)              :: dsdy         ! d(entropy)/(electron fraction)

REAL(KIND=double)              :: denom        ! denominator in Cramer's Law
REAL(KIND=double)              :: drodye_ps    ! d(density)/d(electron fraction)
REAL(KIND=double)              :: drods_pye    ! d(density)/d(entropy)
REAL(KIND=double)              :: dsdr         ! d(entropy)/d(radius)
REAL(KIND=double)              :: dyedr        ! d(electron fraction)/d(radius)

REAL(KIND=double)              :: s_im1        ! enropy at radial zone i-1
REAL(KIND=double)              :: s_ip1        ! enropy at radial zone i+1

REAL(KIND=double)              :: w2BV         ! square of the Brunt-Vaisala frequency
REAL(KIND=double)              :: wBV          ! Brunt-Vaisala frequency (positive if unstable) [s^{-1}]
REAL(KIND=double)              :: twBV         ! Brunt-Vaisala growth time (positive if unstable) [s]

!-----------------------------------------------------------------------
! Shifted radial index
!-----------------------------------------------------------------------
WRITE (*,3002)
 3002 FORMAT (' Entering Brunt_Vaisala_ye')
j                  = i + 1

!-----------------------------------------------------------------------
!  Load EOS scalars
!-----------------------------------------------------------------------

dpdd               = aesvd(j  ,1,ij_ray,ik_ray)
dpdt               = aesvt(j  ,1,ij_ray,ik_ray)
dpdy               = aesvy(j  ,1,ij_ray,ik_ray)
s_im1              = aesv (j-1,3,ij_ray,ik_ray)
s_ip1              = aesv (j+1,3,ij_ray,ik_ray)
dsdd               = aesvd(j  ,3,ij_ray,ik_ray)
dsdt               = aesvt(j  ,3,ij_ray,ik_ray)
dsdy               = aesvy(j  ,3,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Compute d(rho)/d(s) and d(rho)/d(ye)
!-----------------------------------------------------------------------
!
!                       s  drho |
!        drods_pye  :  --- ---- |
!                      rho  ds  |ye,p
!
!                      ye  drho |
!        drodye_ps  :  --- ---- |
!                      rho dye  |s,p
!
!-----------------------------------------------------------------------

denom              = dsdd * dpdt - dsdt * dpdd
drodye_ps          = ( dpdy * dsdt - dpdt * dsdy )/denom
drods_pye          = dpdt/denom

!-----------------------------------------------------------------------
!  Compute dsdr and dyedr
!-----------------------------------------------------------------------

dsdr               = ( s_ip1                   - s_im1                 )/( x_c(i+1) - x_c(i-1) )
dyedr              = ( ye_c(i+1,ij_ray,ik_ray) - ye_c(i-1,ij_ray,ik_ray) )/( x_c(i+1) - x_c(i-1) )

!-----------------------------------------------------------------------
!  Compute the Brunt-Vaisala frequency and time scale
!-----------------------------------------------------------------------

w2BV               = - ( grav_x_c(i,ij_ray,ik_ray)/rho_c(i,ij_ray,ik_ray) ) &
&                  * ( drods_pye * dsdr + drodye_ps * dyedr )
wBV                = SIGN( DSQRT( DABS(w2BV) ), w2BV )
twBV               = 1.d+0/wBV

RETURN
END SUBROUTINE Brunt_Vaisala_ye


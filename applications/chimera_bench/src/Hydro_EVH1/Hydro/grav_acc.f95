SUBROUTINE grav_acc( imin, imax, nx, jmin, jmax, ny, x_e, x_c, mass_a_ray, &
& grav_pot, grav_x )
!-----------------------------------------------------------------------
!
!    File:         grav_acc
!    Module:       grav_acc
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To calculate the gravitational acceleration at the zone centers.
!
!    Input arguments:
!
!  imin         : minimum x-array index
!  imax         : maximum x-array index
!  nx           : x-array extent
!  jmin         : minimum yx-array index
!  jmax         : maximum y-array index
!  nx           : y-array extent
!  x_e          : radial face of zone (cm)
!  x_c          : radial midpoint of zone (cm)
!  mass_a_ray   : mass/angular ray (g)
!
!    Output arguments:
!
!  grav_pot     : gravitational potential
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  e_advct_module, edit_module, evh1_global, evh1_sweep
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY : zero, half, third
USE physcnst_module, ONLY : g
     
IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin          ! minimum unshifted radial zone index
INTEGER, INTENT(in)              :: imax          ! maximum unshifted radial zone index
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: jmin          ! minimum y-zone index
INTEGER, INTENT(in)              :: jmax          ! maximum y-zone index
INTEGER, INTENT(in)              :: ny            ! y-array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)   :: x_e            ! radial face of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: x_c            ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mass_a_ray     ! mass/(angular ray) (g)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,ny) :: grav_pot ! gravitational potential (ergs g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ny) :: grav_x   ! gravitational acceleration (ergs g^{-1} s^{-1})

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: i             ! radial do index
INTEGER                          :: j             ! angular do index

REAL(KIND=double), DIMENSION(nx) :: mass_c        ! enclused mass to center of zone
REAL(KIND=double), DIMENSION(nx) :: chi_f         ! correction factor
REAL(KIND=double), DIMENSION(nx) :: xcen          ! zone center

!-----------------------------------------------------------------------
!  Gravitational potential
!-----------------------------------------------------------------------

mass_c                   = zero
mass_c(imin)             = half * mass_a_ray(imin)

DO i = imin+1,imax
  mass_c(i)              = mass_c(i-1) + half * ( mass_a_ray(i-1) + mass_a_ray(i) )
END DO ! i

DO j = jmin,jmax
  grav_pot(imin:imax,j)  = -g * mass_c(imin:imax)/x_c(imin:imax)
END DO ! j

DO i = imin,imax
  xcen(i)                = 0.5d0 * ( x_e(i+1)**3 + x_e(i)**3 )
  xcen(i)                = DSIGN( DABS(xcen(i))**third, xcen(i) )
END DO

chi_f(imin:imax)         = half * ( x_e(imin:imax)**2 + x_e(imin+1:imax+1)**2 ) &
&                        * half * ( x_e(imin:imax) + x_e(imin+1:imax+1) )/( xcen(imin:imax)**3 )
WRITE (*,3001) (i,chi_f(i),i=imin,imax)
 3001 FORMAT (' i=',i4,' chi_f(i)=',es11.3)

DO j = jmin,jmax
  grav_x(imin:imax,j)    = -g * chi_f * mass_c(imin:imax)/( xcen(imin:imax) * xcen(imin:imax) )
END DO ! j

RETURN
END SUBROUTINE grav_acc

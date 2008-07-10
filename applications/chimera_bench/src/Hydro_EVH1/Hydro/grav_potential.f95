SUBROUTINE grav_potential( imin, imax, nx, jmin, jmax, ny, kmin, kmax, &
& nz, x_e, mass_a_ray, grav_pot )
!-----------------------------------------------------------------------
!
!    File:         grav_potential
!    Module:       grav_potential
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/05
!
!    Purpose:
!      To calculate the gravitational potential at the zone centers.
!
!    Input arguments:
!
!  imin         : minimum x-array index
!  imax         : maximum x-array index
!  nx           : x-array extent
!  jmin         : minimum yx-array index
!  jmax         : maximum y-array index
!  ny           : y-array extent
!  kmin         : minimum zx-array index
!  kmax         : maximum z-array index
!  nz           : z-array extent
!  x_e          : radial midpoint of zone (cm)
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
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY : zero, half
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
INTEGER, INTENT(in)              :: kmin          ! minimum z-zone index
INTEGER, INTENT(in)              :: kmax          ! maximum z-zone index
INTEGER, INTENT(in)              :: nz            ! z-array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: x_e            ! radial midpoint of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)     :: mass_a_ray     ! mass/(angular ray) (g)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,ny,nz) :: grav_pot ! gravitational potential (ergs g^{-1})

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: i             ! radial do index
INTEGER                          :: j             ! angular do index
INTEGER                          :: k             ! azimuthal do index

REAL(KIND=double), DIMENSION(nx) :: mass_e        ! enclused mass to edge of zone
REAL(KIND=double), DIMENSION(nx) :: grav_pot_e    ! mean gravitational potential at zone edge (ergs g^{-1})
REAL(KIND=double), DIMENSION(nx) :: grav_pot_c    ! mean gravitational potential at zone edge (ergs g^{-1})

!-----------------------------------------------------------------------
!  Gravitational potential
!-----------------------------------------------------------------------

mass_e                      = zero

DO i = imin,imax
  mass_e(i+1)               = mass_e(i) + mass_a_ray(i)
END DO ! i

grav_pot_e(imin)            = zero
grav_pot_e(imin+1:imax+1)   = - g * mass_e(imin+1:imax+1)/x_e(imin+1:imax+1)

grav_pot_c(imin:imax)       = half * ( grav_pot_e(imin:imax) + grav_pot_e(imin+1:imax+1) )

DO k = kmin,kmax
  DO j = jmin,jmax
    grav_pot(imin:imax,j,k) = grav_pot_c(imin:imax)
  END DO ! j
END DO ! k

RETURN
END SUBROUTINE grav_potential

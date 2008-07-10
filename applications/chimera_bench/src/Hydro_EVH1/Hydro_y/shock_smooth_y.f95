SUBROUTINE shock_smooth_y( jmin, jmax, ji_ray, jk_ray, ny, j_ray_dim, &
& ik_ray_dim, flat_x_y, rho, j_shock )
!-----------------------------------------------------------------------
!
!    File:         shock_smooth_y
!    Module:       shock_smooth_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To mark zones lying along a shock for additional diffusion
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jmin           : lower y-array index
!  jmax           : upper y-array index
!  ny             : y-array extent
!  ji_ray         : x (radial) index of a specific y (angular) ray
!  jk_ray         : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim      : the number of radial zones on a processor after swapping with y
!  ik_ray_dim     : the number of z-zones on a processor before swapping with z
!  flat_x_y       : variables indicating the presence of radial shocks
!  rho            : densities (g cm^{-3}) after Lagrangian update
!
!    Output arguments:
!  j_shock        : zones marked for additional diffusion
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, parallel_module
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)            :: jmin             ! minimum y-array index
INTEGER, INTENT(in)            :: jmax             ! maximum y-array index
INTEGER, INTENT(in)            :: ny               ! y-array extent

INTEGER, INTENT(in)            :: ji_ray           ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)            :: jk_ray           ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)            :: j_ray_dim        ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)            :: ik_ray_dim       ! number of radial zones on a processor before swapping with z

REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim) :: flat_x_y    ! variables indicating the presence of radial shocks
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim) :: rho         ! density (cm^{-3})

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, DIMENSION(ny)         :: j_shock          ! zones marked for added y-diffusion

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                        :: j                ! y-zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

j_shock            = 0

!-----------------------------------------------------------------------
!  Return if ny < 5
!-----------------------------------------------------------------------

IF ( ny < 5 ) RETURN

!-----------------------------------------------------------------------
!  Check for shock aligned along y-axis, mark zones if odd-even pattern
!   is present
!-----------------------------------------------------------------------

DO j = jmin,jmax-4
  IF ( flat_x_y(j  ,ji_ray,jk_ray) > zero  .and. &
&      flat_x_y(j+1,ji_ray,jk_ray) > zero  .and. &
&      flat_x_y(j+2,ji_ray,jk_ray) > zero  .and. &
&      flat_x_y(j+3,ji_ray,jk_ray) > zero  .and. &
&      flat_x_y(j+4,ji_ray,jk_ray) > zero )        THEN
    IF ( ( rho(j+1,ji_ray,jk_ray) > rho(j  ,ji_ray,jk_ray) + 1.d-6  .and.     &
&          rho(j+2,ji_ray,jk_ray) < rho(j+1,ji_ray,jk_ray) - 1.d-6  .and.     &
&          rho(j+3,ji_ray,jk_ray) > rho(j+2,ji_ray,jk_ray) + 1.d-6  .and.     &
&          rho(j+4,ji_ray,jk_ray) < rho(j+3,ji_ray,jk_ray) - 1.d-6          ) &
&       .or.                                                  &
&        ( rho(j+1,ji_ray,jk_ray) < rho(j  ,ji_ray,jk_ray) - 1.d-6  .and.     &
&          rho(j+2,ji_ray,jk_ray) > rho(j+1,ji_ray,jk_ray) + 1.d-6  .and.     &
&          rho(j+3,ji_ray,jk_ray) < rho(j+2,ji_ray,jk_ray) - 1.d-6  .and.     &
&          rho(j+4,ji_ray,jk_ray) > rho(j+3,ji_ray,jk_ray) + 1.d-6          ) ) THEN
      j_shock(j  ) = 1
      j_shock(j+1) = 1
      j_shock(j+2) = 1
      j_shock(j+3) = 1
    END IF ! rho
  END IF ! flat_x_y
END DO ! j

RETURN
END SUBROUTINE shock_smooth_y

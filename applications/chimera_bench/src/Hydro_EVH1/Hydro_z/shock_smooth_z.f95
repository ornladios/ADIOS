SUBROUTINE shock_smooth_z( kmin, kmax, ki_ray, kj_ray, nz, ij_ray_dim, &
& k_ray_dim, flat_x_z, rho, k_shock )
!-----------------------------------------------------------------------
!
!    File:         shock_smooth_z
!    Module:       shock_smooth_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/12/07
!
!    Purpose:
!      To mark zones lying along a shock for additional diffusion
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  kmin           : lower z-array index
!  kmax           : upper z-array index
!  nz             : z-array extent
!  ki_ray         : x (radial) index of a specific z (azimuthal) ray
!  kj_ray         : z (azimuthal) index of a specific z (azimuthal) ray
!  ij_ray_dim     : the number of radial zones on a processor before swapping with y
!  k_ray_dim      : the number of z-zones on a processor after swapping with z
!  flat_x_z       : variables indicating the presence of radial shocks
!  rho            : densities (g cm^{-3}) after Lagrangian update
!
!    Output arguments:
!  k_shock        : zones marked for additional diffusion
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, parallel_module
!
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

INTEGER, INTENT(in)            :: kmin             ! minimum y-array index
INTEGER, INTENT(in)            :: kmax             ! maximum y-array index
INTEGER, INTENT(in)            :: nz               ! y-array extent

INTEGER, INTENT(in)            :: ki_ray           ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)            :: kj_ray           ! y (azimuthal) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)            :: ij_ray_dim       ! number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)            :: k_ray_dim        ! number of radial zones on a processor after swapping with z

REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim) :: flat_x_z    ! variables indicating the presence of radial shocks
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim) :: rho         ! density (cm^{-3})

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, DIMENSION(nz)         :: k_shock          ! zones marked for added y-diffusion

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                        :: iz               ! z-zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

k_shock            = 0

!-----------------------------------------------------------------------
!  Return if nz < 5
!-----------------------------------------------------------------------

IF ( nz < 5 ) RETURN

!-----------------------------------------------------------------------
!  Check for shock aligned along y-axis, mark zones if odd-even pattern
!   is present
!-----------------------------------------------------------------------

DO iz = kmin,kmax-4
  IF ( flat_x_z(iz  ,kj_ray,ki_ray) > zero  .and. &
&      flat_x_z(iz+1,kj_ray,ki_ray) > zero  .and. &
&      flat_x_z(iz+2,kj_ray,ki_ray) > zero  .and. &
&      flat_x_z(iz+3,kj_ray,ki_ray) > zero  .and. &
&      flat_x_z(iz+4,kj_ray,ki_ray) > zero )        THEN
    IF ( ( rho(iz+1,kj_ray,ki_ray) > rho(iz  ,kj_ray,ki_ray) + 1.d-6  .and.     &
&          rho(iz+2,kj_ray,ki_ray) < rho(iz+1,kj_ray,ki_ray) - 1.d-6  .and.     &
&          rho(iz+3,kj_ray,ki_ray) > rho(iz+2,kj_ray,ki_ray) + 1.d-6  .and.     &
&          rho(iz+4,kj_ray,ki_ray) < rho(iz+3,kj_ray,ki_ray) - 1.d-6          ) &
&       .or.                                                  &
&        ( rho(iz+1,kj_ray,ki_ray) < rho(iz  ,kj_ray,ki_ray) - 1.d-6  .and.     &
&          rho(iz+2,kj_ray,ki_ray) > rho(iz+1,kj_ray,ki_ray) + 1.d-6  .and.     &
&          rho(iz+3,kj_ray,ki_ray) < rho(iz+2,kj_ray,ki_ray) - 1.d-6  .and.     &
&          rho(iz+4,kj_ray,ki_ray) > rho(iz+3,kj_ray,ki_ray) + 1.d-6          ) ) THEN
      k_shock(iz  ) = 1
      k_shock(iz+1) = 1
      k_shock(iz+2) = 1
      k_shock(iz+3) = 1
    END IF ! rho
  END IF ! flat_x_z
END DO ! z = kmin,kmax-4

RETURN
END SUBROUTINE shock_smooth_z

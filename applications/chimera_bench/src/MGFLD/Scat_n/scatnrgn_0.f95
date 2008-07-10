SUBROUTINE scatnrgn_0( j, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         scatnrgn_0
!    Module:       scatnrgn_0
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/09/03
!
!    Purpose:
!      To set to zero the nucleon elastic scattering functions at the
!       corners of the "unit cell" in rho-T-ye state space.
!
!    Input arguments:
!
!  j        : radial zone number
!  ij_ray   : j-index of a radial ray
!  ik_ray   : k-index of a radial ray
!
!    Include files:
!  array_module
!  nu_energy_grid_module, scat_n_module
!
!-----------------------------------------------------------------------

USE array_module, ONLY : ij_ray_dim

USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE scat_n_module, ONLY : sctn0, sctnb0

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: k             ! incomiong neutrino energy zone index
INTEGER                          :: kp            ! outgoing neutrino energy zone index
INTEGER                          :: ida           ! density cube index
INTEGER                          :: ita           ! temperature cube index
INTEGER                          :: iya           ! electron fraction cube index

INTEGER                          :: i_ray         ! f(ij_ray,ik_ray)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Combine indices to uniquely index the radial ray
!-----------------------------------------------------------------------

i_ray              = ij_ray_dim * ( ik_ray - 1 ) + ij_ray

!-----------------------------------------------------------------------
!  Set rates to zero
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  DO kp = 1,k
    DO ida = 1,2
      DO ita = 1,2
        DO iya = 1,2
          sctn0 (j,k,kp,ida,ita,iya,i_ray) = -100.e0
          sctn0 (j,kp,k,ida,ita,iya,i_ray) = -100.e0
          sctnb0(j,k,kp,ida,ita,iya,i_ray) = -100.e0
          sctnb0(j,kp,k,ida,ita,iya,i_ray) = -100.e0
        END DO
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE scatnrgn_0

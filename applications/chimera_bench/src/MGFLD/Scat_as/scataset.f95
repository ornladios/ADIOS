SUBROUTINE scataset( j, ij_ray, ik_ray, rho, t, ye )
!-----------------------------------------------------------------------
!
!    File:         scataset
!    Module:       scataset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To set the neutrino-nucleon inelastic scattering rates.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  j             : shifted radial zone index
!  ij_ray        : j-index of a radial ray
!  ik_ray        : k-index of a radial ray
!  rho           : density (cm^{-3})
!  t             : temperature (K)
!  ye            : electron fraction
!
!    Output arguments:
!        none
!
!    Include files:
!      kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: j             ! radial zone index
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)     :: rho           ! density (g cm^{-3})
REAL(KIND=double), INTENT(in)     :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)     :: ye            ! electron fraction

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE scataset

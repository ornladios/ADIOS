SUBROUTINE sctarate( n, jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         sctarate
!    Module:       sctarate
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To evaluate the neutrino-nucleon inelastic scattering rates.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!        none
!
!-----------------------------------------------------------------------

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor index
INTEGER, INTENT(in)               :: jr_min        ! minimum radial index
INTEGER, INTENT(in)               :: jr_max        ! maximum radial index
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE sctarate

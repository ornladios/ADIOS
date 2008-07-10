SUBROUTINE luminosity( jr_min, jr_max, i_ray, i_ray_dim, nx, nnu, lum )
!-----------------------------------------------------------------------
!
!    File:         luminosity
!    Module:       luminosity
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/12/00
!
!    Purpose:
!      To calculate the neutrino luminosities at each radial zone interface.
!
!    Subprograms called:
!  flux            : Computes the neutrino flux as a funciton of j and n
!
!    Input arguments:
!  jr_min          : inner radial zone number
!  jr_max          : outer radial zone number
!  i_ray           : index denoting a specific radial ray
!  i_ray_dim       : number radial rays assigned to a processor
!  nx              : x-array extent
!  nnu             : neutrino flavor array extent
!
!    Output arguments:
!  lum             : luminosity at each zone interface
!
!    Include files:
!  kind_module, array_module, physcnst_module
!  edit_module, mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpi
USE physcnst_module, ONLY : ergfoe, cvel

USE edit_module, ONLY : nprint
USE mdl_cnfg_module, ONLY : r
USE nu_dist_module, ONLY : fluxnu
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: jr_min        ! minimum radial index
INTEGER, INTENT(in)               :: jr_max        ! maximum radial index
INTEGER, INTENT(in)               :: i_ray         ! index denoting a specific radial ray
INTEGER, INTENT(in)               :: i_ray_dim     ! number radial rays assigned to a processor
INTEGER, INTENT(in)               :: nx            ! x-array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,nnu,i_ray_dim) :: lum    ! luminosity at radius r_lum

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: j             ! radial zone index
INTEGER                           :: n             ! neutrino flavor index

REAL(KIND=double)                 :: area_j        ! area at j

!-----------------------------------------------------------------------
!        Initialize.
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO j = jr_min,jr_max
    lum(j-jr_min+1,n,i_ray) = zero
  END DO
END DO


!-----------------------------------------------------------------------
!        Return if no neutrinos.
!-----------------------------------------------------------------------

IF (nnugpmx == 0) RETURN

!-----------------------------------------------------------------------
!
!                \\\\\ NEUTRINO LUMINOSITIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Compute the luminosities as a function of j and n.
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  CALL flux( jr_min, jr_max, n )
  DO j = jr_min,jr_max
    area_j                  = frpi * r(j) * r(j)
    lum(j-jr_min+1,n,i_ray) = area_j * fluxnu(j,n) * ergfoe
  END DO ! j
END DO ! n

RETURN
END SUBROUTINE luminosity

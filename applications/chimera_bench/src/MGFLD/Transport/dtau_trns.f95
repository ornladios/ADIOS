SUBROUTINE dtau_trns( jr_min, jr_max, dtnph_trans )
!-----------------------------------------------------------------------
!
!    File:         dtau_trns
!    Module:       dtau_trns
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/09/00
!
!    Purpose:
!      To compute the timestep for neutrino absorption, emission,
!       and transport.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  jr_min          : inner radial zone number
!  jr_max          : outer radial zone number
!  dtnph_trans     : coordinate time step for neutrino source and transport
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  agr_nu(j)       : lapse function for mass zone j
!
!    Output arguments (common):
!
!  dtauj_nutrns(j) : proper time step for neutrino absorption, emission,
!  dtau_nutrns(j)  : proper time step for neutrino absorption, emission,
!                    and transport in radial zone j - 1/2
!  cdt(j)          : !*a*dtau_nutrns(j)
!  cdtinv(j)       : 1/!*a*dtau_nutrns(j)
!
!    Include files:
!  numerical_module, physcnst_module
!  nu_dist_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : half, one
USE physcnst_module, ONLY : cvel

USE nu_dist_module, ONLY : agr_nu
USE t_cntrl_module, ONLY : dtauj_nutrns, dtau_nutrns, cdt, cdtinv

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)              :: jr_min          ! minimum radial zone index
INTEGER, INTENT(IN)              :: jr_max          ! maximum radial zone index

REAL(KIND=double), INTENT(in)    :: dtnph_trans     ! source and transport time step

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index

!-----------------------------------------------------------------------
!       Compute cdt * agr(j) and its inverse, and the time step for
!        for emission, absorption, and transport.
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  dtauj_nutrns(j)  = dtnph_trans * agr_nu(j)
  dtau_nutrns(j)   = dtnph_trans * ( half * ( agr_nu(j) + agr_nu(j-1) ) )
  cdt(j)           = cvel * dtau_nutrns(j)
  cdtinv(j)        = one/cdt(j)
END DO

RETURN
END SUBROUTINE dtau_trns

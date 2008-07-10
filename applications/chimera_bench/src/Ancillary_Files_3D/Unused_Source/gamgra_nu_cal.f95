SUBROUTINE gamgra_nu_cal( jr_min, jr_max )
!-----------------------------------------------------------------------
!
!    File:         gamgra_nu_cal
!    Module:       gamgra_nu_cal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/10/00
!
!    Purpose:
!      To compute the GR gamma to be used in neutrino transport.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  jr_min      : inner zone for which calculation agr is to be made
!  jr_max      : outer zone for which calculation agr is to be made
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  gamgr(j)    : inverse radial metric function
!  ireltrns    : GR transport switch
!
!    Output arguments (common):
!
! gamgra_nu(j) : inverse radial metric function at time m to be
!                 in the neutrino transport modules
!
!    Include files:
!  numerical_module
!  edit_module, mdl_cnfg.cmn, nu_dist.cmn, prb_cntl.cmn!
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : one

USE edit_module, ONLY : nprint
USE mdl_cnfg_module, ONLY : gamgr
USE nu_dist_module, ONLY : gamgra_nu
USE prb_cntl_module, ONLY : ireltrns

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(IN)              :: jr_max        ! maximum radial zone index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  201 FORMAT (' ireltrns =',i5,' not supported by gamgra_nu_cal')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

GR: SELECT CASE (ireltrns)

!-----------------------------------------------------------------------
!
!               \\\\\ NONRELATIVISTIC TRANSPORT /////
!
!-----------------------------------------------------------------------

 CASE(0,1)

  DO j = jr_min-1,jr_max+1
    gamgra_nu(j)     = one
  END DO

  RETURN

!-----------------------------------------------------------------------
!
!                \\\\\ RELATIVISTIC TRANSPORT /////
!
!-----------------------------------------------------------------------

 CASE(2)

  DO j = jr_min-1,jr_max+1
    gamgra_nu(j)     = gamgr(j)
  END DO
  
  RETURN

!-----------------------------------------------------------------------
!  Stop if irelhy /= 0 and irelhy /= 1
!-----------------------------------------------------------------------

 CASE DEFAULT GR

  WRITE (nprint, 201) ireltrns
  STOP

END SELECT GR

END SUBROUTINE gamgra_nu_cal
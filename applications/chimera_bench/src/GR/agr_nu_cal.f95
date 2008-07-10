SUBROUTINE agr_nu_cal( jr_min, jr_max )
!-----------------------------------------------------------------------
!
!    File:         agr_nu_cal
!    Module:       agr_nu_cal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/15/03
!
!    Purpose:
!      To compute the GR lapse function at time m to be used in the
!       neutrino transport modules.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min      : inner zone for which calculation agr is to be made
!  jr_max      : outer zone for which calculation agr is to be made
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  agr(j)      : lapse function
!  ireltrns    : GR transport switch
!
!    Output arguments (common):
!  agr_nu(j)   : lapse function at time m to be used in the neutrino transport modules
!
!    Include files:
!  numerical_module
!  edit_module, mdl_cnfg.cmn, nu_dist.cmn, prb_cntl.cmn
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY: one

USE edit_module, ONLY : nprint
USE mdl_cnfg_module, ONLY: agr, agrh
USE nu_dist_module, ONLY: agr_nu, agrjmh_nu
USE prb_cntl_module, ONLY: ireltrns

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(IN)              :: jr_max        ! maximum radial zone index

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  201 FORMAT (' ireltrns =',i5,' not supported by agr_nu_cal')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

GR: SELECT CASE (ireltrns)

!-----------------------------------------------------------------------
!
!               \\\\\ NONRELATIVISTIC TRANSPORT /////
!
!-----------------------------------------------------------------------

 CASE(0)

  agr_nu   (jr_min-1:jr_max+1)     = one
  agrjmh_nu(jr_min-1:jr_max+1)     = one
  
  RETURN

!-----------------------------------------------------------------------
!
!                \\\\\ RELATIVISTIC TRANSPORT /////
!
!-----------------------------------------------------------------------

 CASE(1)

  agr_nu   (jr_min-1:jr_max+1)     = agr(jr_min-1:jr_max+1)
  agrjmh_nu(jr_min-1:jr_max)       = agrh(jr_min-1:jr_max)
  
  RETURN

!-----------------------------------------------------------------------
!  Stop if irelhy /= 0 and irelhy /= 1
!-----------------------------------------------------------------------

 CASE DEFAULT GR

  WRITE (nprint, 201) ireltrns
  STOP

END SELECT GR

END SUBROUTINE agr_nu_cal

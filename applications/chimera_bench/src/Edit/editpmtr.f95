SUBROUTINE editpmtr
!-----------------------------------------------------------------------
!
!    File:         editpmtr
!    Module:       editpmtr
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/2/00
!
!    Purpose:
!      To modify edit control parameters during a run.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!      none
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!
!  nedx(i)    : x-edit counter for data set i.
!  intedx(i)  : number of cycles between x-edits of data set i.
!
!    Include files:
!  cycle_module, edit_module
!
!-----------------------------------------------------------------------

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : intdsc, intedy

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE editpmtr

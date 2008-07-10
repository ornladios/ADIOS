SUBROUTINE model_check( i_radhyd_data, i_model_data )
!-----------------------------------------------------------------------
!
!    File:         model_check
!    Module:       model_check
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/04
!
!    Purpose:
!      Check data for consistency.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  i_radhyd_data   : integer array of radhyd keys
!  i_model_data    : integer array of initial model data
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE edit_module, ONLY : nprint, nlog

IMPLICIT none
SAVE


!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in), DIMENSION(50)           :: i_radhyd_data ! integer array of edit keys
INTEGER, INTENT(in), DIMENSION(2)            :: i_model_data  ! integer array of initial model data

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                  :: imin       ! minimum x-array index
INTEGER                                  :: imax       ! maximum x-array index
INTEGER                                  :: imin_mdl   ! minimum x-array index given by model
INTEGER                                  :: imax_mdl   ! maximum x-array index given by model

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' ***WARNING*** imin_mdl =',i6,' /= 1')
 1003 FORMAT (' imin =',i6,' < imin_mdl =',i6,'(declared lower radial extent is less than data')
 1005 FORMAT (' imin =',i6,' > imin_mdl =',i6,'(declared lower radial extent is greater than data')
 1007 FORMAT (' imax =',i6,' < imax_mdl =',i6,'(declared upper radial extent is less than data')
 1009 FORMAT (' imax =',i6,' > imax_mdl =',i6,'(declared upper radial extent is greater than data')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

imin                = i_radhyd_data(13)
imax                = i_radhyd_data(14)
imin_mdl            = i_model_data(1)
imax_mdl            = i_model_data(2)

!-----------------------------------------------------------------------
!
!               \\\\\ CHECK CONSISTENCY OF DATA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Test prescribed extents against model data extents
!-----------------------------------------------------------------------

IF ( imin_mdl > 1 ) THEN
  WRITE (nlog,1001) imin_mdl
  WRITE (nprint,1001) imin_mdl
END IF ! jr_max + 1 > nx

IF ( imin_mdl > imin ) THEN
  WRITE (nlog,1003) imin, imin_mdl
  WRITE (nprint,1003) imin, imin_mdl
  STOP
END IF ! jr_max + 1 > nx

IF ( imin_mdl < imin ) THEN
  WRITE (nlog,1005) imin, imin_mdl
  WRITE (nprint,1005) imin, imin_mdl
  STOP
END IF ! jr_max + 1 > nx

IF ( imax_mdl > imax ) THEN
  WRITE (nlog,1007) imax, imax_mdl
  WRITE (nprint,1007) imax, imax_mdl
  STOP
END IF ! jr_max + 1 > nx

IF ( imax_mdl < imax ) THEN
  WRITE (nlog,1009) imax, imax_mdl
  WRITE (nprint,1009) imax, imax_mdl
  STOP
END IF ! jr_max + 1 > nx


RETURN
END SUBROUTINE model_check

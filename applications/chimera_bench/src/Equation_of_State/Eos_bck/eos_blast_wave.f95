SUBROUTINE eos( i_ray )

!***********************************************************************

USE kind_module, ONLY : double
USE eos_bck_module, ONLY : jshel, dbck, tbck, yebck, erad, prad, enu, pnu, sneu, &
& dtran, etot, eh, ee, ed, egy0, ptotbck, ph, pe, pd, stotbck, sh, se, sd, dedt
USE eos_snc_x_module, ONLY : nse

USE numerical_module, ONLY : zero
USE physcnst_module, ONLY : dmnp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: i_ray         ! index denoting a specific radial ray

REAL(KIND=double)                :: srad          ! entropy (nucleon^{-1})


!*************************************************** radiation *********

erad            = zero
prad            = zero
srad            = zero
enu             = erad
pnu             = prad
sneu            = srad

CALL eos0

!************************************************** electrons *********
!*************************************************no neutrinos ********

CALL lectron

!**************************************************** nucleons *********

CALL eosnuc_x

!************************************************ get totals ***********

etot            = ed + ee
ptotbck         = pd + pe
stotbck         = sd + se

!***************************************** approx for energy inversion *

dedt            = ed/tbck

RETURN
END

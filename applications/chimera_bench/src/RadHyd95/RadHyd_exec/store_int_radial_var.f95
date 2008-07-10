 SUBROUTINE store_int_radial_var
!-----------------------------------------------------------------------
!
!    File:         store_int_radial_var
!    Module:       store_int_radial_var
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To store coordinate and state variables before the Lagrangian hydro step.
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
!  numerical_module
!  revh1_global, incrmnt_module, adial_ray_module
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE evh1_global, ONLY : svel
USE incrmnt_module, ONLY : dtmpmn, dtmpnn, dye, dyecnvt
USE radial_ray_module, ONLY : rho_c, t_c, ye_c, ei_c, u_c, v_c, w_c, &
& rho_ci, t_ci, ye_ci, ei_ci, u_ci, v_ci, w_ci
     
IMPLICIT none

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!         \\\\\ STORE INITIAL VARIABLES FOR THE X-SWEEP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  State variables
!-----------------------------------------------------------------------

rho_ci               = rho_c
t_ci                 = t_c
ye_ci                = ye_c
ei_ci                = ei_c
u_ci                 = u_c
v_ci                 = v_c
w_ci                 = w_c

!-----------------------------------------------------------------------
!  Increment variables
!-----------------------------------------------------------------------

dtmpmn               = zero
dtmpnn               = zero
dye                  = zero
dyecnvt              = zero

!-----------------------------------------------------------------------
!  Courant variavle
!-----------------------------------------------------------------------

svel                 = zero

RETURN
END SUBROUTINE store_int_radial_var

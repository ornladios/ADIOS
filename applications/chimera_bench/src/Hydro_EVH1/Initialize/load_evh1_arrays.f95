SUBROUTINE load_evh1_arrays
!-----------------------------------------------------------------------
!
!    File:         load_evh1_arrays
!    Module:       load_evh1_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!        To load evh1 arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!    Output arguments:
!        none
!
!    Include files:
!      numerical_module
!      radial_ray_module, evh1_global, evh1_bound
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE radial_ray_module, ONLY : ndimp=>ndim, ngeomxp=>ngeomx, ngeomyp=>ngeomy, &
& ngeomzp=>ngeomz, nleftxp=> nleftx, nleftyp=>nlefty, nleftzp=>nleftz, &
& nrightxp=>nrightx, nrightyp=>nrighty, nrightzp=>nrightz


USE radial_ray_module, ONLY : u_l=>u_bcl, v_l=>v_bcl, w_l=>w_bcl, r_l=>r_bcl, &
& p_l=>p_bcl, ei_l=>ei_bcl, ye_l=>ye_bcl, temp_l=>temp_bcl, gc_l=>gc_bcl, ge_l=>ge_bcl, &
& psi0_l=>psi0_bcl, comp_l=>comp_bcl, u_r=>u_bcr, v_r=>v_bcr, w_r=>w_bcr, r_r=>r_bcr, &
& p_r=>p_bcr, ei_r=>ei_bcr, ye_r=>ye_bcr, temp_r=>temp_bcr, gc_r=>gc_bcr, ge_r=>ge_bcr, &
& psi0_r=>psi0_bcr, comp_r=>comp_bcr

USE evh1_global, ONLY : ndim, ngeomx, ngeomy, ngeomz, nleftx, nlefty, nleftz, nrightx, &
& nrighty, nrightz

USE evh1_bound, ONLY : u_bcl, v_bcl, w_bcl, r_bcl, p_bcl, ei_bcl, ye_bcl, temp_bcl, &
& gc_bcl, ge_bcl, psi0_bcl, comp_bcl, u_bcr, v_bcr, w_bcr, r_bcr, p_bcr, ei_bcr, &
& ye_bcr, temp_bcr, gc_bcr, ge_bcr, psi0_bcr, comp_bcr

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Geometry
!-----------------------------------------------------------------------

ndim                   = ndimp
ngeomx                 = ngeomxp
ngeomy                 = ngeomyp
ngeomz                 = ngeomzp
nleftx                 = nleftxp
nlefty                 = nleftyp
nleftz                 = nleftzp
nrightx                = nrightxp
nrighty                = nrightyp
nrightz                = nrightzp

!-----------------------------------------------------------------------
!  Left boundary values
!-----------------------------------------------------------------------

u_bcl                  = u_l
v_bcl                  = v_l
w_bcl                  = w_l
r_bcl                  = r_l
p_bcl                  = p_l
ei_bcl                 = ei_l
ye_bcl                 = ye_l
temp_bcl               = temp_l
gc_bcl                 = gc_l
ge_bcl                 = ge_l
psi0_bcl               = psi0_l
comp_bcl               = comp_l

!-----------------------------------------------------------------------
!  Right boundary values
!-----------------------------------------------------------------------

u_bcr                  = u_r
v_bcr                  = v_r
w_bcr                  = w_r
r_bcr                  = r_r
p_bcr                  = p_r
ei_bcr                 = ei_r
ye_bcr                 = ye_r
temp_bcr               = temp_r
gc_bcr                 = gc_r
ge_bcr                 = ge_r
psi0_bcr               = psi0_r
comp_bcr               = comp_r


RETURN
END SUBROUTINE load_evh1_arrays

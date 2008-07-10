SUBROUTINE model_write( ndump, nx, nnu )
!-----------------------------------------------------------------------
!
!    File:         model_write
!    Module:       model_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/25/05
!
!    Purpose:
!      To dump the model configuration.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ndump       : unit number to dump the data
!  nx          : x-array extent
!  nnu         : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  edit_modulee, eos_snc_x_module, nu_dist_module, nu_energy_grid_module,
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : nlog, nu_r, nu_rt, nu_rho, nu_rhot, psi0dat, psi1dat
USE eos_snc_x_module, ONLY : duesrc
USE nu_dist_module, ONLY : dnurad, unukrad, unujrad, nnukrad, nnujrad, &
& e_rad, unurad, elec_rad, nnurad
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE radial_ray_module, ONLY : imin, imax, ncycle, nprint, nse_c, rho_c, &
& t_c, ye_c, u_c, v_c, w_c, x_ef, dx_cf, psi0_c, xn_c, a_nuc_rep_c, &
& z_nuc_rep_c, be_nuc_rep_c, uburn_c, e_nu_c_bar, f_nu_e_bar

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ndump           ! unit number to write restart file
INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Formats: Document the dump
!-----------------------------------------------------------------------

 1001 FORMAT (' ***Model dump written at cycle      ',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!
!               \\\\\ WRITE THE MODEL RESTART FILES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Radial index bounds for MGFLD shifted radial arrays
!-----------------------------------------------------------------------

WRITE (ndump) imin
WRITE (ndump) imax

!-----------------------------------------------------------------------
!  Independent thermodynamic variables
!-----------------------------------------------------------------------

WRITE (ndump) rho_c
WRITE (ndump) t_c
WRITE (ndump) ye_c

!-----------------------------------------------------------------------
!  Independent mechanical variables
!-----------------------------------------------------------------------

WRITE (ndump) u_c
WRITE (ndump) v_c
WRITE (ndump) w_c
WRITE (ndump) dx_cf
WRITE (ndump) x_ef

!-----------------------------------------------------------------------
!  Independent radiation variables and bookkeeping arrays
!-----------------------------------------------------------------------

WRITE (ndump) psi0_c
WRITE (ndump) dnurad
WRITE (ndump) unukrad
WRITE (ndump) unujrad
WRITE (ndump) e_rad
WRITE (ndump) unurad
WRITE (ndump) nnukrad
WRITE (ndump) nnujrad
WRITE (ndump) elec_rad
WRITE (ndump) nnurad
WRITE (ndump) e_nu_c_bar
WRITE (ndump) f_nu_e_bar

!-----------------------------------------------------------------------
!  Net number of neutrinos radiated from density rho_nurad and radius
!   r_nurad
!-----------------------------------------------------------------------

WRITE (ndump) nu_r
WRITE (ndump) nu_rt
WRITE (ndump) nu_rho
WRITE (ndump) nu_rhot

!-----------------------------------------------------------------------
!  Time integrated psi0 and psi1
!-----------------------------------------------------------------------

WRITE (ndump) psi0dat
WRITE (ndump) psi1dat

!-----------------------------------------------------------------------
!  nse - non-bse flag
!-----------------------------------------------------------------------

WRITE (ndump) nse_c

!-----------------------------------------------------------------------
!  Nuclear abundances
!-----------------------------------------------------------------------

WRITE (ndump) xn_c

!-----------------------------------------------------------------------
!  Auxiliary heavy nucleus
!-----------------------------------------------------------------------

WRITE (ndump) a_nuc_rep_c
WRITE (ndump) z_nuc_rep_c
WRITE (ndump) be_nuc_rep_c

!-----------------------------------------------------------------------
!  Nuclear energy released
!-----------------------------------------------------------------------

WRITE (ndump) uburn_c

!-----------------------------------------------------------------------
!  Energy offsets
!-----------------------------------------------------------------------

WRITE (ndump) duesrc

!-----------------------------------------------------------------------
!
!                    \\\\\ RECORD THE DUMP /////
!
!-----------------------------------------------------------------------

WRITE (nlog,1001) ncycle,ndump

RETURN
END SUBROUTINE model_write

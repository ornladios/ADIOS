SUBROUTINE dimension_arrays( nx, ny, nz, nez, nnu, nnc, max_12, n_proc, &
& ij_ray_dim, ik_ray_dim, j_ray_dim, k_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_arrays
!    Module:       dimension_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/04
!
!    Purpose:
!      To dimension the array dimensions.
!
!    Subprograms called:
!  dimension_radhyd_arrays
!  dimension_hydro_arrays
!  dimension_mgfld_arrays
!  dimension_nucbrn_arrays
!  dimension_e_advct_arrays
!
!    Input arguments:
!  nx         : x (radial) array dimension
!  ny         : y (angular) array dimension
!  nz         : z (azimuthal) array dimension
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!  nnc        : composition array dimension
!  max_12     : max(nx,ny,nz)+12
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  j_ray_dim  : number of radial zones on a processor after swapping with y
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!
!    Output arguments:
!        none
!
!    Include files:
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nx            ! x (radial array) dimension
INTEGER, INTENT(in)               :: ny            ! y (angular) array dimension
INTEGER, INTENT(in)               :: nz            ! z (azimuthal) array dimension
INTEGER, INTENT(in)               :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array dimension
INTEGER, INTENT(in)               :: nnc           ! composition array dimension
INTEGER, INTENT(in)               :: max_12        ! max(nx,ny,nz)+12
INTEGER, INTENT(in)               :: n_proc        ! number of processors assigned to run
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)               :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)               :: k_ray_dim     ! number of radial zones on a processor after swapping with z

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Calling dimension_radhyd_arrays')
  103 FORMAT (' Radhyd arrays have been dimensioned')
  105 FORMAT (' Calling dimension_hydro_arrays')
  107 FORMAT (' Hydro arrays have been dimensioned')
  109 FORMAT (' Calling dimension_tov_pot_arrays')
  111 FORMAT (' GR arrays have been dimensioned')
  113 FORMAT (' Calling dimension_mgfld_arrays')
  115 FORMAT (' MGFLD arrays have been dimensioned')
  117 FORMAT (' Calling dimension_edit_arrays')
  119 FORMAT (' Edit arrays have been dimensioned')
  121 FORMAT (' Calling dimension_eos_arrays')
  123 FORMAT (' Equation of state arrays have been dimensioned')
  125 FORMAT (' Calling dimension_nucbrn_arrays')
  127 FORMAT (' Nuclear arrays have been dimensioned')
  129 FORMAT (' Calling dimension_e_advct_arrays')
  131 FORMAT (' Neutrino energy advection arrays have been dimensioned')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Set radhyd array dimensions.
!
!        Calls
!         dimension_radial_var_arrays
!         dimension_angular_ray_arrays
!         dimension_prb_cntl_arrays
!         dimension_t_cntrl_arrays
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,101)
CALL dimension_radhyd_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, j_ray_dim, &
& k_ray_dim, nez, nnu, nnc )
IF ( myid == 0 ) WRITE (nlog,103)

!-----------------------------------------------------------------------
!        Set hydro array dimensions.
!
!        Calls
!         dimension_boundary_arrays
!         dimension_convect_arrays
!         dimension_mgfld_remap_arrays
!         dimension_shock_arrays
!         dimension_evh1_sweep_arrays
!         dimension_evh1_zone_arrays
!         dimension_evh1_bound_arrays
!         dimension_mdl_cnfg_y_arrays
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,105)
CALL dimension_hydro_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, j_ray_dim, &
& k_ray_dim, nez, nnu, nnc, max_12 )
IF ( myid == 0 ) WRITE (nlog,107)

!-----------------------------------------------------------------------
!        Set GR array dimensions.
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,109)
CALL dimension_tov_pot_arrays( nx )
IF ( myid == 0 ) WRITE (nlog,111)

!-----------------------------------------------------------------------
!        Set mgfld array dimensions.
!
!        Calls
!         dimension_abem_arrays
!         dimension_abem_y_arrays
!         dimension_brem_arrays
!         dimension_incrmnt_arrays
!         dimension_mdl_cnfg_arrays
!         dimension_nu_dist_arrays
!         dimension_nu_energy_grid_arrays
!         dimension_pair_arrays
!         dimension_scat_a_arrays
!         dimension_scat_e_arrays
!         dimension_scat_i_arrays
!         dimension_scat_n_arrays
!         dimension_scat_nA_arrays
!         dimension_scat_nn_arrays
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,113)
CALL dimension_mgfld_arrays( nx, ny, nz, nez, nnu, nnc, ij_ray_dim, &
& ik_ray_dim, j_ray_dim, k_ray_dim )
IF ( myid == 0 ) WRITE (nlog,115)

!-----------------------------------------------------------------------
!        Set edit array dimensions.
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,117)
CALL dimension_edit_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, nez, nnu )
IF ( myid == 0 ) WRITE (nlog,119)

!-----------------------------------------------------------------------
!        Set equation of statge array dimensions.
!
!        Calls
!         dimension_eos_snc_x_arrays
!         dimension_eos_snc_y_arrays
!         dimension_eos_snc_z_arrays
!         dimension_eos_bck_arrays
!         dimension_eos_ls_arrays
!         dimension_eos_drv_arrays
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,121)
CALL dimension_eos_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, j_ray_dim, &
& k_ray_dim, nnc )
IF ( myid == 0 ) WRITE (nlog,123)

!-----------------------------------------------------------------------
!        Set nuclear array dimensions.
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,125)
CALL dimension_nucbrn_arrays( nx, nnc, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,127)

!-----------------------------------------------------------------------
!        Set e_advection array dimensions.
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,129)
CALL dimension_e_advct_arrays( nx, ny, nz, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,131)

RETURN
END SUBROUTINE dimension_arrays

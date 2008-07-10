SUBROUTINE dimension_mgfld_arrays( nx, ny, nz, nez, nnu, nnc, ij_ray_dim, &
& ik_ray_dim, j_ray_dim, k_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_mgfld_arrays
!    Module:       dimension_mgfld_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To allocate array dimensions.
!
!    Subprograms called:
!  dimension_abem_arrays           : dimensions x - emmision and absorption arrays
!  dimension_abem_y_arrays         : dimensions y - emmision and absorption arrays
!  dimension_brem_arrays           : dimensions bremsstrahlung arrays
!  dimension_mdl_cnfg_arrays       : dimensions model x - configuration arrays
!  dimension_mdl_cnfg_y_arrays     : dimensions model y - configuration arrays
!  dimension_mdl_cnfg_z_arrays     : dimensions model z - configuration arrays
!  dimension_nu_dist_arrays        : dimensions neutrino arrays
!  dimension_nu_energy_grid_arrays : dimensions the neutrino energy grid arrays
!  dimension_pair_arrays           : dimensions electron-positron pair annihilation arrays
!  dimension_pair_A_arrays         : dimensions the nuclear deexcitation pair arrays
!  dimension_scat_a_arrays         : dimension the inelastic neutrino-nuclerus scattering arrays (Haxton)
!  dimension_scat_e_arrays         : dimension the neutrino-electron scattering arrays
!  dimension_scat_i_arrays         : dimension the isoenergetic scattering arrays
!  dimension_scat_n_arrays         : dimension the neutrino-nucleon scattering arrays
!  dimension_scat_n_arrays         : dimension the neutrino-nucleon elastic scattering arrays
!  dimension_scat_nA_arrays        : dimension the neutrino-nucleus inelastic scattering arrays
!  dimension_scat_nn_arrays        : dimension the neutrino-nucleon inelastic scattering arrays
!
!    Input arguments:
!  nx            : x (radial) array dimension
!  ny            : y (angular) array dimension
!  nz            : z (asimuthal) array dimension
!  nez           : energy array extent
!  nnu           : neutrino flavor array dimension
!  nnc           : composition array dimension
!  ij_ray_dim    : number of y-zones on a processor before swapping with y
!  ik_ray_dim    : number of z-zones on a processor before swapping with z
!  j_ray_dim     : number of radial zones on a processor after swapping with y
!  k_ray_dim     : number of radial zones on a processor after swapping with z
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

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nx         ! x (radial) array extent
INTEGER, INTENT(in)               :: ny         ! y (angular) array extent
INTEGER, INTENT(in)               :: nz         ! y (asimuthal) array extent
INTEGER, INTENT(in)               :: nez        ! neutrino energy array extent
INTEGER, INTENT(in)               :: nnu        ! neutrino flavor array extent
INTEGER, INTENT(in)               :: nnc        ! composition array dimension
INTEGER, INTENT(in)               :: ij_ray_dim ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)               :: ik_ray_dim ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)               :: j_ray_dim  ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)               :: k_ray_dim  ! number of radial zones on a processor after swapping with z

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Calling dimension_abem_arrays')
  103 FORMAT (' Calling dimension_abem_y_arrays')
  105 FORMAT (' Calling dimension_abem_z_arrays')
  107 FORMAT (' Calling dimension_brem_arrays')
  109 FORMAT (' Calling dimension_incrmnt_arrays')
  111 FORMAT (' Calling dimension_mdl_cnfg_arrays')
  113 FORMAT (' Calling dimension_mdl_cnfg_y_arrays')
  115 FORMAT (' Calling dimension_mdl_cnfg_z_arrays')
  117 FORMAT (' Calling dimension_nu_dist_arrays')
  119 FORMAT (' Calling dimension_nu_energy_grid_arrays')
  121 FORMAT (' Calling dimension_pair_arrays')
  123 FORMAT (' Calling dimension_pair_A_arrays')
  125 FORMAT (' Calling dimension_scat_a_arrays')
  127 FORMAT (' Calling dimension_scat_e_arrays')
  129 FORMAT (' Calling dimension_scat_i_arrays')
  131 FORMAT (' Calling dimension_scat_n_arrays')
  133 FORMAT (' Calling dimension_scat_nA_arrays')
  135 FORMAT (' Calling dimension_scat_nn_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,101)
CALL dimension_abem_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,103)
CALL dimension_abem_y_arrays( ny, nez, nnu, j_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,105)
CALL dimension_abem_z_arrays( nz, nez, nnu, ij_ray_dim, k_ray_dim )
IF ( myid == 0 ) WRITE (nlog,107)
CALL dimension_brem_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,109)
CALL dimension_incrmnt_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,111)
CALL dimension_mdl_cnfg_arrays( nx )
IF ( myid == 0 ) WRITE (nlog,113)
CALL dimension_mdl_cnfg_y_arrays( ny )
IF ( myid == 0 ) WRITE (nlog,115)
CALL dimension_mdl_cnfg_z_arrays( nz )
IF ( myid == 0 ) WRITE (nlog,117)
CALL dimension_nu_dist_arrays( nx, ny, nz, nez, nnu, nnc, ij_ray_dim, &
& ik_ray_dim, j_ray_dim, k_ray_dim )
IF ( myid == 0 ) WRITE (nlog,119)
CALL dimension_nu_energy_grid_arrays( nez, nnu )
IF ( myid == 0 ) WRITE (nlog,121)
CALL dimension_pair_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,123)
!CALL dimension_pair_A_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,125)
CALL dimension_scat_a_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,127)
CALL dimension_scat_e_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,129)
CALL dimension_scat_i_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,131)
CALL dimension_scat_n_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,133)
CALL dimension_scat_nA_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
IF ( myid == 0 ) WRITE (nlog,135)
CALL dimension_scat_nn_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )

RETURN
END SUBROUTINE dimension_mgfld_arrays

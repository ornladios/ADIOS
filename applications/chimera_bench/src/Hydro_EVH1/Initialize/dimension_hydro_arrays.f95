SUBROUTINE dimension_hydro_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, nez, nnu, nnc, max_12 )
!-----------------------------------------------------------------------
!
!    File:         dimension_hydro_arrays
!    Module:       dimension_hydro_arrays
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
!  dimension_boundary_arrays    : dimensions boundary_module arrays
!  dimension_convect_arrays     : dimensions convect_module arrays
!  dimension_evh1_bound_arrays  : dimensions evh1_bound arrays
!  dimension_evh1_sweep_arrays  : dimensions evh1_sweep arrays
!  dimension_evh1_zone_arrays   : dimensions evh1_zone arrays
!  dimension_mgfld_remap_arrays : dimensions mgfld_remap_module arrays
!  dimension_shock_arrays       : dimensions shock_module arrays
!
!    Input arguments:
!  nx         : x (radial) array extent
!  ny         : y (angular) array extent
!  nz         : z (azimuthal) array extent
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  j_ray_dim  : the number of radial zones on a processor after swapping
!                with y
!  k_ray_dim  : the number of radial zones on a processor after swapping
!                with z
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  nnc        : composition array extent
!  max_12     : max(nx,ny,nz)+12
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

INTEGER, INTENT(in)               :: nx            ! x (radial) array extent
INTEGER, INTENT(in)               :: ny            ! y (angular) array extent
INTEGER, INTENT(in)               :: nz            ! z (azimuthal) array extent
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)               :: j_ray_dim     ! the number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)               :: k_ray_dim     ! the number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)               :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array dimension
INTEGER, INTENT(in)               :: nnc           ! composition array dimension
INTEGER, INTENT(in)               :: max_12        ! max(nx,ny,nz)+12

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Calling dimension_boundary_arrays')
  103 FORMAT (' Calling dimension_convect_arrays')
  105 FORMAT (' Calling dimension_mgfld_remap_arrays')
  107 FORMAT (' Calling dimension_shock_arrays')
  109 FORMAT (' Calling dimension_evh1_sweep_arrays')
  111 FORMAT (' Calling dimension_evh1_zone_arrays')
  113 FORMAT (' Calling dimension_evh1_bound_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,101)
CALL dimension_boundary_arrays( nx )
IF ( myid == 0 ) WRITE (nlog,103)
CALL dimension_convect_arrays( nx )
IF ( myid == 0 ) WRITE (nlog,105)
CALL dimension_mgfld_remap_arrays( max_12, nez, nnu, nnc )
IF ( myid == 0 ) WRITE (nlog,107)
CALL dimension_shock_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim )
IF ( myid == 0 ) WRITE (nlog,109)
CALL dimension_evh1_sweep_arrays( max_12 )
IF ( myid == 0 ) WRITE (nlog,111)
CALL dimension_evh1_zone_arrays( nx, ny, nz )
IF ( myid == 0 ) WRITE (nlog,113)
CALL dimension_evh1_bound_arrays( nnc )

RETURN
END SUBROUTINE dimension_hydro_arrays

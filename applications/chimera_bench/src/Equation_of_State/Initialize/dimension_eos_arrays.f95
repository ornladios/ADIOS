SUBROUTINE dimension_eos_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, nnc )
!-----------------------------------------------------------------------
!
!    File:         dimension_eos_arrays
!    Module:       dimension_eos_arrays
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
!  dimension_eos_snc_x_arrays : allocate dimensions to eos_snc_x_module arrays 
!  dimension_eos_snc_y_arrays : allocate dimensions to eos_snc_y_module arrays 
!  dimension_eos_bck_arrays   : allocate dimensions to eos_bck_module arrays 
!  dimension_eos_ls_arrays    : allocate dimensions to eos_ls_module arrays
!  dimension_eos_drv_arrays   ! allocate dimensions to eos_drv_module arrays
!
!    Input arguments:
!  nx         : x (radial) array extent
!  ny         : y (angular) array extent
!  nz         : z (azimuthal) array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  j_ray_dim  : number of angular rays assigned to a processor
!  nnc        : composition array dimension
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

INTEGER, INTENT(in)               :: nx            ! x-array extent
INTEGER, INTENT(in)               :: ny            ! y-array extent
INTEGER, INTENT(in)               :: nz            ! z-array extent
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)               :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)               :: k_ray_dim     ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)               :: nnc           ! composition array dimension

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Calling dimension_eos_snc_x_arrays')
  103 FORMAT (' Calling dimension_eos_snc_y_arrays')
  105 FORMAT (' Calling dimension_eos_snc_z_arrays')
  107 FORMAT (' Calling dimension_eos_bck_arrays')
  109 FORMAT (' Calling dimension_eos_ls_arrays')
  111 FORMAT (' Calling dimension_eos_drv_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,101)
CALL dimension_eos_snc_x_arrays( nx, ij_ray_dim, ik_ray_dim, nnc )
IF ( myid == 0 ) WRITE (nlog,103)
CALL dimension_eos_snc_y_arrays( ny, j_ray_dim, ik_ray_dim, nnc )
IF ( myid == 0 ) WRITE (nlog,105)
CALL dimension_eos_snc_z_arrays( nz, ij_ray_dim, k_ray_dim, nnc )
IF ( myid == 0 ) WRITE (nlog,107)
CALL dimension_eos_bck_arrays( nx )
IF ( myid == 0 ) WRITE (nlog,109)
CALL dimension_eos_ls_arrays( nx )
IF ( myid == 0 ) WRITE (nlog,111)
CALL dimension_eos_drv_arrays( nx, ij_ray_dim, ik_ray_dim )

RETURN
END SUBROUTINE dimension_eos_arrays

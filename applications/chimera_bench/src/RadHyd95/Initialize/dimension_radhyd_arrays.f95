SUBROUTINE dimension_radhyd_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         dimension_radhyd_arrays
!    Module:       dimension_radhyd_arrays
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
!      dimension_radial_var_arrays
!      dimension_radhyd_ray_arrays
!      dimension_angular_ray_arrays
!      dimension_prb_cntl_arrays
!      dimension_t_cntrl_arrays
!
!    Input arguments:
!  nx         : x (radial) array extent
!  ny         : y angular array extent
!  nz         : z azimuthal array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  j_ray_dim  : number of radial zones on a processor after swapping with y
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
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
INTEGER, INTENT(in)               :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array dimension
INTEGER, INTENT(in)               :: nnc           ! composition array dimension

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Calling dimension_radial_ray_arrays')
  103 FORMAT (' Calling dimension_angular_ray_arrays')
  105 FORMAT (' Calling dimension_azimuthal_ray_arrays')
  107 FORMAT (' Calling dimension_prb_cntl_arrays')
  109 FORMAT (' Calling dimension_t_cntrl_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( myid == 0 ) WRITE (nlog,101)
CALL dimension_radial_ray_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, nez, nnu, nnc )
IF ( myid == 0 ) WRITE (nlog,103)
CALL dimension_angular_ray_arrays( ny, j_ray_dim, ik_ray_dim, nez, nnu, &
& nnc )
IF ( myid == 0 ) WRITE (nlog,105)
CALL dimension_azimuthal_ray_arrays( nz, ij_ray_dim, k_ray_dim, nez, &
& nnu, nnc )
IF ( myid == 0 ) WRITE (nlog,107)
CALL dimension_prb_cntl_arrays( nnu )
IF ( myid == 0 ) WRITE (nlog,109)
CALL dimension_t_cntrl_arrays( nx, nnu )

RETURN
END SUBROUTINE dimension_radhyd_arrays

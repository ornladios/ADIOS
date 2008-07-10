SUBROUTINE dimension_shock_arrays( nx, ny, nz, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_shock_arrays
!    Module:       dimension_shock_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the shock arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x (radial) array extent
!  ny         : y (angular) array extent
!  nz         : z (azimuthal) array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  j_ray_dim  : number of radial zones on a processor after swapping with y
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, parallel_module, shock_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE shock_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x (radial array) dimension
INTEGER, INTENT(in)              :: ny            ! y (angular array) dimension
INTEGER, INTENT(in)              :: nz            ! z (azimuthal array) dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: k_ray_dim     ! number of radial zones on a processor after swapping with z

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Shock arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_shock_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!              \\\\\ ALLOCATE SHOCK_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (lshock(nx+ny+nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'lshock    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (pq_x(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pq_x      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pqr_x(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pqr_x     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pqy_x(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pqy_x     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pqz_x(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pqz_x     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (q0_x(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'q0_x      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (pq_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pq_y      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pqr_y(ny,j_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pqr_y     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (q0_y(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'q0_y      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (pq_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pq_y      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pqr_z(nz,ij_ray_dim,k_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pqr_y     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (q0_z(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'q0_y      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (j_shk_radial_p(ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'j_shk__p  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (j_shk_radial_all_p(ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'j_shk_all '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (nqpmax(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nqpmax    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (qpmax(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qpmax     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (qrpmax(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qrpmax    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!              \\\\\ INITIALIZE SHOCK_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

ipq                       = 0

pq_x                      = zero
pqr_x                     = zero
pqy_x                     = zero
q0_x                      = -1.d0
pqcrit                    = zero

pq_y                      = zero
pqr_y                     = zero
q0_y                      = -1.d0

pq_z                      = zero
pqr_z                     = zero
q0_z                      = -1.d0

jshockmn                  = 0
jshockmx                  = 0
j_shk_radial_p            = 0
j_shk_radial_all_p        = 0

nqpmax                    = 0
qpmax                     = zero
qrpmax                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_shock_arrays

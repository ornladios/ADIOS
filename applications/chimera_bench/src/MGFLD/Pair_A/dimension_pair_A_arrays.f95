SUBROUTINE dimension_pair_A_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_pair_A_arrays
!    Module:       dimension_pair_A_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the nuclear pair annihilation and production arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x (radial) array dimension
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, pair_A_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE pair_A_module
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Nuclear pair annihilation and production arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_pair_A_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!              \\\\\ ALLOCATE PAIR_A_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrp_A(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrp_A    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrp_A(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrp_A    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrp_A(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrp_A    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Nuclear pair annihilation and production functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (pair_A_a(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pair_A_a  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pair_A_p(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pair_A_p  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated Nuclear pair annihilation and production annihilation
!   functions
!-----------------------------------------------------------------------

ALLOCATE (paf_A(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paf_A     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafd_A(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafd_A    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (paft_A(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paft_A    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafy_A(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafy_A    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (pafp0_A(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafp0_A   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafp1_A(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafp1_A   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ INITIALIZE PAIR_A_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrp_A                    = 0
itrp_A                    = 0
iyrp_A                    = 0

!-----------------------------------------------------------------------
!  Nuclear pair annihilation and production functions at cube corners
!-----------------------------------------------------------------------

pair_A_a                  = -100.d0
pair_A_p                  = -100.d0

!-----------------------------------------------------------------------
!  Interpolated Nuclear pair annihilation and production annihilation
!   functions
!-----------------------------------------------------------------------

paf_A                     = zero
pafd_A                    = zero
paft_A                    = zero
pafy_A                    = zero
pafp0_A                   = zero
pafp1_A                   = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_pair_A_arrays

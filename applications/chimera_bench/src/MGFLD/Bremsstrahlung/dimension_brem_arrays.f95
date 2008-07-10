SUBROUTINE dimension_brem_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_brem_arrays
!    Module:       dimension_brem_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the bremsstrahlung arrays.
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
!  brem_module, edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE brem_module
USE edit_module, ONLY : nlog
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

  101 FORMAT (' Radhyd time and time step control arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_brem_arrays')

!-----------------------------------------------------------------------
!
!              \\\\\ ALLOCATE BREM_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrb(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrb      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrb(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrb      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrb(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrb      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Bremsstrahlung pair annihilation functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (brema0(nx,nez,nez,2,2,2,ij_ray_dim * ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'brema0    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated bremsstrahlung pair annihilation functions
!-----------------------------------------------------------------------

ALLOCATE (baf(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'baf       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (bafd(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'bafd      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (baft(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'baft      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (bafy(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'bafy      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (bafp0(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'bafp0     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (bafp1(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'bafp1     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!              \\\\\ INITIALIZE BREM_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrb                      = 0
itrb                      = 0
iyrb                      = 0

!-----------------------------------------------------------------------
!  Bremsstrahlung pair annihilation functions at cube corners
!-----------------------------------------------------------------------

brema0                    = -100d0

!-----------------------------------------------------------------------
!  Interpolated bremsstrahlung pair annihilation functions
!-----------------------------------------------------------------------

baf                       = zero
bafd                      = zero
baft                      = zero
bafy                      = zero
bafp0                     = zero
bafp1                     = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_brem_arrays

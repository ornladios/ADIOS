SUBROUTINE dimension_pair_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_pair_arrays
!    Module:       dimension_pair_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the pair production arrays.
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
!  edit_module, pair_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE pair_module
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

  101 FORMAT (' Electron-positron pair annihilation arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_pair_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!               \\\\\ ALLOCATE PAIR_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrpp(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrpp     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrpp(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrpp     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrpp(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrpp     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (paira0i(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paira0i   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (paira0ii(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paira0ii  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated bremsstrahlung pair annihilation functions
!-----------------------------------------------------------------------

ALLOCATE (paf(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paf       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafd(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafd      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (paft(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paft      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafy(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafy      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (pafp0(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafp0     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafp1(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafp1     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE PAIR_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrpp                     = 0
itrpp                     = 0
iyrpp                     = 0

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation functions at cube corners
!-----------------------------------------------------------------------

paira0i                   = -100.d0
paira0ii                  = -100.d0

!-----------------------------------------------------------------------
!  Interpolated bremsstrahlung pair annihilation functions
!-----------------------------------------------------------------------

paf                       = zero
pafd                      = zero
paft                      = zero
pafy                      = zero
pafp0                     = zero
pafp1                     = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_pair_arrays

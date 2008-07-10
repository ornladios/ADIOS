SUBROUTINE dimension_scat_nn_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_scat_nn_arrays
!    Module:       dimension_scat_nn_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the neutrino-nucleon inelastic scattering arrays.
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
!  edit_module, parallel_module, scat_nn_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE scat_nn_module

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

  101 FORMAT (' Neutrino-nucleon inelastic scattering arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_scat_nn_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE SCAT_NN_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrnn(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrnn     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrnn(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrnn     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrnn(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrnn     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino-nucleon inelastic scattering functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (sctnn0(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sctnn0    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated Neutrino-nucleon inelastic scattering functions
!-----------------------------------------------------------------------

ALLOCATE (scnnf(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnnf     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnnfd(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnnfd    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnnft(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnnft    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnnfy(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnnfy    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (scnnfp0(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnnfp0   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnnfp1(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnnfp1   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE SCAT_NN_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrnn                     = 0
itrnn                     = 0
iyrnn                     = 0

!-----------------------------------------------------------------------
!  Neutrino-nucleon inelastic scattering functions at cube corners
!-----------------------------------------------------------------------

sctnn0                    = -100.d0

!-----------------------------------------------------------------------
!  Interpolated Neutrino-nucleon inelastic scattering functions
!-----------------------------------------------------------------------

scnnf                     = zero
scnnfd                    = zero
scnnft                    = zero
scnnfy                    = zero

scnnfp0                   = zero
scnnfp1                   = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_scat_nn_arrays

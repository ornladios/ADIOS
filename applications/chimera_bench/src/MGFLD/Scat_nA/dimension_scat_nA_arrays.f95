SUBROUTINE dimension_scat_nA_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim)
!-----------------------------------------------------------------------
!
!    File:         dimension_scat_nA_arrays
!    Module:       dimension_scat_nA_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the neutrino-nucleon inelastic
!       scattering arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x (radial) array dimension
!  nez       : neutrino energy array dimension
!  nnu       : neutrino flavor array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, parallel_module, scat_nA_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE scat_nA_module

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

  101 FORMAT (' Neutrino-nucleus inelastic scattering arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_scat_nA_arrays')

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

ALLOCATE (idrnA(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrnA     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrnA(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrnA     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrnA(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrnA     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino-nucleon inelastic scattering functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (sctnA0(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sctnA0    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated Neutrino-nucleon inelastic scattering functions
!-----------------------------------------------------------------------

ALLOCATE (scnAf(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnAf     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnAfd(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnAfd    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnAft(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnAft    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnAfy(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnAfy    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (scnAfp0(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnAfp0   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnAfp1(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnAfp1   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE SCAT_NN_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrnA                     = 0
itrnA                     = 0
iyrnA                     = 0

!-----------------------------------------------------------------------
!  Neutrino-nucleon inelastic scattering functions at cube corners
!-----------------------------------------------------------------------

sctnA0                    = -100.d0

!-----------------------------------------------------------------------
!  Interpolated Neutrino-nucleon inelastic scattering functions
!-----------------------------------------------------------------------

scnAf                     = zero
scnAfd                    = zero
scnAft                    = zero
scnAfy                    = zero
scnAfp0                   = zero
scnAfp1                   = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_scat_nA_arrays

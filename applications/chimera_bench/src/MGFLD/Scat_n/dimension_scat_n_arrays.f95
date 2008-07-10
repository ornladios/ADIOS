SUBROUTINE dimension_scat_n_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_scat_n_arrays
!    Module:       dimension_scat_n_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the neutrino-nucleon elastic scattering arrays.
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
!  edit_module, parallel_module, scat_n_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE scat_n_module

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

  101 FORMAT (' Neutrino-nucleon elastic scattering arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_scat_n_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE SCAT_N_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrn(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrn      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrn(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrn      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrn(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrn      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (sctn0(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sctn0     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Antieutrino-nucleon elastic scattering functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (sctnb0(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sctnb0    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated neutrino-nucleon inelastic scattering functions
!-----------------------------------------------------------------------

ALLOCATE (scnf(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnf      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnfd(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnfd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnft(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnft     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnfy(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnfy     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (scnfp0(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnfp0    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnfp1(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnfp1    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE SCAT_N_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrn                      = 0
itrn                      = 0
iyrn                      = 0

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering functions at cube corners
!-----------------------------------------------------------------------

sctn0                     = -100.d0

!-----------------------------------------------------------------------
!  Antieutrino-nucleon elastic scattering functions at cube corners
!-----------------------------------------------------------------------

sctnb0                    = -100.d0

!-----------------------------------------------------------------------
!  Interpolated neutrino-nucleon inelastic scattering functions
!-----------------------------------------------------------------------

scnf                      = zero
scnfd                     = zero
scnft                     = zero
scnfy                     = zero
scnfp0                    = zero
scnfp1                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_scat_n_arrays

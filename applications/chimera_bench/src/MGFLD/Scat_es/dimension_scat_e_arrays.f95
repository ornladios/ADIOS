SUBROUTINE dimension_scat_e_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_scat_e_arrays
!    Module:       dimension_scat_e_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the neutrino-electron scattering arrays.
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
!  edit_module, parallel_module, scat_e_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE scat_e_module

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

  101 FORMAT (' Neutrino-electron elastic scattering arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_scat_e_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE SCAT_E_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrse(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrse     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrse(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrse     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrse(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrse     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Neutrino-electron scattering functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (scte0i(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scte0i    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scte0ii(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scte0ii   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated neutrino-electron scattering functions
!-----------------------------------------------------------------------

ALLOCATE (scef(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scef      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scefd(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scefd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (sceft(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sceft     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scefy(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scefy     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (scefp0(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scefp0    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scefp1(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scefp1    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE SCAT_E_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrse                     = 0
itrse                     = 0
iyrse                     = 0

!-----------------------------------------------------------------------
!  Neutrino-electron scattering functions at cube corners
!-----------------------------------------------------------------------

scte0i                    = -100.d0
scte0ii                   = -100.d0

!-----------------------------------------------------------------------
!  Interpolated neutrino-electron scattering functions
!-----------------------------------------------------------------------

scef                      = zero
scefd                     = zero
sceft                     = zero
scefy                     = zero
scefp0                    = zero
scefp1                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_scat_e_arrays

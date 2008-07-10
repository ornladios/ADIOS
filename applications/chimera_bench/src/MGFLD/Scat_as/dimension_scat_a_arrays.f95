SUBROUTINE dimension_scat_a_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_scat_a_arrays
!    Module:       dimension_scat_a_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the Haxton rate arrays.
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
!  edit_module, parallel_module, scat_a_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE scat_a_module

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

  101 FORMAT (' The Haxton neutrino-nucleus inelastic scattering arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_scat_a_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE SCAT_A_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrsa(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrsa     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrsa(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrsa     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrsa(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrsa     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  n-neutrino-nucleus non-isoenergetic scattering data
!-----------------------------------------------------------------------

ALLOCATE (rncnu0(nez,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rncnu0    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rncnb0(nez,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rncnb0    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  n-neutrino-nucleus inelastic scattering functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (scta0_nue(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scta0_nue '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scta0_nueb(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scta0_nueb'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scta0_nux(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scta0_nux '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scta0_nuxb(nx,nez,nez,2,2,2,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scta0_nuxb'; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated n-neutrino-nucleus inelastic scattering functions
!-----------------------------------------------------------------------

ALLOCATE (scaf(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scaf      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scafd(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scafd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scaft(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scaft     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scafy(6,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scafy     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (scafp0(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scafp0    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scafp1(6,nez,nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scafp1    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE SCAT_A_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrsa                     = 0
itrsa                     = 0
iyrsa                     = 0

!-----------------------------------------------------------------------
!  n-neutrino-nucleus non-isoenergetic scattering data
!-----------------------------------------------------------------------

rncnu0                    = zero
rncnb0                    = zero

!-----------------------------------------------------------------------
!  n-neutrino-nucleus inelastic scattering functions at cube corners
!-----------------------------------------------------------------------

scta0_nue                 = zero
scta0_nueb                = zero
scta0_nux                 = zero
scta0_nuxb                = zero

!-----------------------------------------------------------------------
!  Interpolated n-neutrino-nucleus inelastic scattering functions
!-----------------------------------------------------------------------

scaf                      = zero
scafd                     = zero
scaft                     = zero
scafy                     = zero
scafp0                    = zero
scafp1                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_scat_a_arrays

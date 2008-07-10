SUBROUTINE dimension_scat_i_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_scat_i_arrays
!    Module:       dimension_scat_i_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the isoenergetic scattering arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x (radial) array dimension
!  nez        : neutrino energy array dimension
!  nez        : neutrino energy array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, parallel_module, scat_i_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE scat_i_module

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

  101 FORMAT (' Neutrino isoenergetic scattering arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_scat_i_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE SCAT_I_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrsi(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrsi     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrsi(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrsi     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrsi(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrsi     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Isoenergetic scattering functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (cohsct(nx,nez,2,2,2,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cohsct    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cohbsct(nx,nez,2,2,2,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cohbsct   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated isoenergetic scattering inverse mean free paths
!-----------------------------------------------------------------------

ALLOCATE (scti(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scti      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (sctid(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sctid     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (sctit(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sctit     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (sctiy(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sctiy     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE SCAT_I_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrsi                     = 0
itrsi                     = 0
iyrsi                     = 0

!-----------------------------------------------------------------------
!  Isoenergetic neutrino scattering functions at cube corners
!-----------------------------------------------------------------------

cohsct                    = -100.d0

!-----------------------------------------------------------------------
!  Isoenergetic antineutrino scattering functions at cube corners
!-----------------------------------------------------------------------

cohbsct                   = -100.d0

!-----------------------------------------------------------------------
!  Interpolated isoenergetic neutrino scattering inverse mean free paths
!-----------------------------------------------------------------------

scti                      = zero
sctid                     = zero
sctit                     = zero
sctiy                     = zero

!-----------------------------------------------------------------------
!  Isoenergetic neutrino scattering inverse mean free paths
!-----------------------------------------------------------------------

rmdnps0                   = zero
rmdnns0                   = zero
rmdnhes0                  = zero
rmdnhs0                   = zero
rmdnps1                   = zero
rmdnns1                   = zero
rmdnhes1                  = zero
rmdnhs1                   = zero

!-----------------------------------------------------------------------
!  Isoenergetic antineutrino scattering inverse mean free paths
!-----------------------------------------------------------------------

rmdnbps0                  = zero
rmdnbns0                  = zero
rmdnbps1                  = zero
rmdnbns1                  = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_scat_i_arrays

SUBROUTINE dimension_abem_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_abem_arrays
!    Module:       dimension_abem_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the aben arrays.
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
!  abem_module, edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE abem_module
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

  101 FORMAT (' Emission and absorption arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_abem_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!              \\\\\ ALLOCATE ABEM_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices.
!-----------------------------------------------------------------------

ALLOCATE (idrae(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrae     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrae(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrae     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrae(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrae     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Haxton's emission and absporption rates
!-----------------------------------------------------------------------

ALLOCATE (rncem(nez,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rncem     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rncep(nez,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rncep     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Totsl emission and absporption inverse mean free paths at cube
!   corners
!-----------------------------------------------------------------------

ALLOCATE (em(nx,nez,nnu,2,2,2,ij_ray_dim * ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'em        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ab(nx,nez,nnu,2,2,2,ij_ray_dim * ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ab        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Emission inverse mean free paths
!-----------------------------------------------------------------------

ALLOCATE (emis(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'emis      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (emisd(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'emisd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (emist(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'emist     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (emisy(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'emisy     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Absorption inverse mean free paths
!-----------------------------------------------------------------------

ALLOCATE (absor(nx,nez,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'absor     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (absord(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'absord    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (absort(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'absort    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (absory(nx,nez,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'absory    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!              \\\\\ INITIALIZE ABEM_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrae                     = 0
itrae                     = 0
iyrae                     = 0

!-----------------------------------------------------------------------
!  Haxton's emission and absporption rates
!-----------------------------------------------------------------------

rncem                     = zero
rncep                     = zero

!-----------------------------------------------------------------------
!  Totsl emission and absporption inverse mean free paths at cube
!   corners
!-----------------------------------------------------------------------

em                        = zero
ab                        = zero

!-----------------------------------------------------------------------
!  Emission inverse mean free paths
!-----------------------------------------------------------------------

emis                      = zero
emisd                     = zero
emist                     = zero
emisy                     = zero

!-----------------------------------------------------------------------
!  Absorption inverse mean free paths
!-----------------------------------------------------------------------

absor                     = zero
absord                    = zero
absort                    = zero
absory                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_abem_arrays

SUBROUTINE dimension_incrmnt_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_incrmnt_arrays
!    Module:       dimension_incrmnt_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the increment arrays.
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
!  edit_module, incrmnt_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE incrmnt_module
USE parallel_module, ONLY : myid

IMPLICIT none

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

CHARACTER (len=11)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Increment arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a11,' in dimension_incrmnt_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE INCRMNT_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature increments
!-----------------------------------------------------------------------

ALLOCATE (dtmpmn(nx,10,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtmpmn     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dtmpnn(nx,10,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtmpnn     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Internal energy increments
!-----------------------------------------------------------------------

ALLOCATE (denergy_int(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'denergy_int'; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Electron fraction increments
!-----------------------------------------------------------------------

ALLOCATE (dye(nx,nnu,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dyecnvt(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dyecnvt    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE INCRMNT_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature increments
!-----------------------------------------------------------------------

dtmpmn                    = zero
dtmpnn                    = zero

!-----------------------------------------------------------------------
!  Internal energy increments
!-----------------------------------------------------------------------

denergy_int               = zero

!-----------------------------------------------------------------------
!  Electron fraction increments
!-----------------------------------------------------------------------

dye                       = zero
dyecnvt                   = zero

WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_incrmnt_arrays

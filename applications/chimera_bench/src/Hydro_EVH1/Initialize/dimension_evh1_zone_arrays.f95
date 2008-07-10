SUBROUTINE dimension_evh1_zone_arrays(nx,ny,nz)
!-----------------------------------------------------------------------
!
!    File:         dimension_evh1_zone_arrays
!    Module:       dimension_evh1_zone_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to certain of the evh1 arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x (radial) array dimension
!            : y (angular array) dimension
!  nz        : z dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  edit_module, evh1_zone, parallel_module
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : nlog
USE evh1_zone
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nx            ! x (radial array) dimension
INTEGER, INTENT(in)               :: ny            ! y (angular array) dimension
INTEGER, INTENT(in)               :: nz            ! z dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' EVH1 zone arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_evh1_zone_arrays')

!-----------------------------------------------------------------------
!        Allocate evh1_zone_module arrays
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                     \\\\\ STATE VARIABLES /////
!
!-----------------------------------------------------------------------

ALLOCATE (zro(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zro       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zei(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zei       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zegrav(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zegrav    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zeg(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zeg       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zux(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zux       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zuy(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zuy       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zuz(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zuz       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zte(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zte       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zye(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zye       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (znu_str(nx,ny,nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'znu_str   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                     \\\\\ COORDINATE GRID /////
!
!-----------------------------------------------------------------------

ALLOCATE (zxa(nx+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zxa       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zya(ny+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zya       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zza(nz+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zza       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zdx(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zdx       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zdy(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zdy       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zdz(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zdz       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zxc(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zxc       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zyc(ny), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zyc       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zzc(nz), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zzc       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!              \\\\\ PPM INTERPOLATION COEFFICIENTS /////
!
!-----------------------------------------------------------------------

ALLOCATE (zparax(10,nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zparax    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zparay(10,ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zparay    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zparaz(10,nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zparaz    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                      \\\\\ INITIALIZATION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  State Variables
!-----------------------------------------------------------------------

zro                       = zero
zei                       = zero
zegrav                    = zero
zeg                       = zero
zux                       = zero
zuy                       = zero
zuz                       = zero
zte                       = zero
zye                       = zero
znu_str                   = zero

!-----------------------------------------------------------------------
!  Coordinate grid
!-----------------------------------------------------------------------

zxa                       = zero
zya                       = zero
zza                       = zero
zdx                       = zero
zdy                       = zero
zdz                       = zero
zxc                       = zero
zyc                       = zero
zzc                       = zero

!-----------------------------------------------------------------------
!  PPM Interpolation coefficients
!-----------------------------------------------------------------------

zparax                    = zero
zparay                    = zero
zparaz                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_evh1_zone_arrays
                                                                                

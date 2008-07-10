SUBROUTINE dimension_eos_drv_arrays( nx, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_eos_drv_arrays
!    Module:       dimension_eos_drv_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/02/08
!
!    Purpose:
!      To allocate dimensions to eos_drv_module arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x-array extent
!  ij_ray_dim,ik_ray_dim  : number of radial rays  assigned to a processor
!  nnc        : composition aray extent
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, eos_drv_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE eos_drv_module
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: ij_ray_dim,ik_ray_dim     ! number of radial rays to put on a processor

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Equation of state arrays for radial rays have been dimensioned in dimension_eos_drv_arrays')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_eos_drv_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!            \\\\\ ALLOCATE EOS_SNC_X_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Interpolated equation of state variables
!-----------------------------------------------------------------------

ALLOCATE (aesv(nx,3,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesv      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvd(nx,3,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvd     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvt(nx,3,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvt     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (aesvy(nx,3,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aesvy     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idr(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itr(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyr(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyr       '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  EOS table
!-----------------------------------------------------------------------

ALLOCATE (estble(3,nx,2,2,2,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'estble    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (escnst(3,nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'escnst    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Density derivative EOS table
!-----------------------------------------------------------------------

ALLOCATE (estbled(3,nx,2,2,2,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'estbled   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (escnstd(3,nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'escnstd   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Temperature derivative EOS table
!-----------------------------------------------------------------------

ALLOCATE (estblet(3,nx,2,2,2,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'estblet   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (escnstt(3,nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'escnstt   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Lepton fraction derivative EOS table
!-----------------------------------------------------------------------

ALLOCATE (estbley(3,nx,2,2,2,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'estbley   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (escnsty(3,nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'escnsty   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Lepton fraction array
!-----------------------------------------------------------------------

ALLOCATE (yl(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yl        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!           \\\\\ INITIALIZE EOS_SNC_X_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Interpolated equation of state variables
!-----------------------------------------------------------------------

aesv                      = zero
aesvd                     = zero
aesvt                     = zero
aesvy                     = zero

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idr                       = 0
itr                       = 0
iyr                       = 0

!-----------------------------------------------------------------------
!  EOS table
!-----------------------------------------------------------------------

estble                    = zero
escnst                    = zero

!-----------------------------------------------------------------------
!  Density deerivative EOS table
!-----------------------------------------------------------------------

estbled                   = zero
escnstd                   = zero

!-----------------------------------------------------------------------
!  Temperature deerivative EOS table
!-----------------------------------------------------------------------

estblet                   = zero
escnstt                   = zero

!-----------------------------------------------------------------------
!  Lepton fraction deerivative EOS table
!-----------------------------------------------------------------------

estbley                   = zero
escnsty                   = zero

!-----------------------------------------------------------------------
!  Lepton fraction array
!-----------------------------------------------------------------------

yl                        = zero

!-----------------------------------------------------------------------
!  Density above which matter and neutrinos are coupled
!-----------------------------------------------------------------------

rho_couple                = 1.d+12

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_eos_drv_arrays

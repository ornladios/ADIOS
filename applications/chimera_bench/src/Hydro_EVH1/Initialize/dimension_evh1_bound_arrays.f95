SUBROUTINE dimension_evh1_bound_arrays(nnc)
!-----------------------------------------------------------------------
!
!    File:         dimension_evh1_bound_arrays
!    Module:       dimension_evh1_bound_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate the dimensions certain of the evh1 arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x (radial) array dimension
!  ny        : y (angular array) dimension
!  nz        : z dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, evh1_bound, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE evh1_bound
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nnc           ! number of species

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                :: var_name

INTEGER                           :: istat         ! allocation status

  101 FORMAT (' evh1_bound_arrays have been initialized')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_evh1_sweep_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Allocate evh1_bound module arrays
!-----------------------------------------------------------------------

ALLOCATE (comp_bcl(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp_bcl  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (comp_bcr(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp_bcr  '; WRITE (nlog,1001) var_name; END IF

!........Initialize evh1_bound module arrays............................
!.......................................................................
!.......................................................................

!........Left boundary initialization...................................
!.......................................................................

u_bcl                     = zero
v_bcl                     = zero
w_bcl                     = zero
r_bcl                     = zero
p_bcl                     = zero
ei_bcl                    = zero
ye_bcl                    = zero
temp_bcl                  = zero
gc_bcl                    = zero
ge_bcl                    = zero
psi0_bcl                  = zero
comp_bcl                  = zero

!........Left boundary initialization...................................
!.......................................................................

u_bcr                     = zero
v_bcr                     = zero
w_bcr                     = zero
r_bcr                     = zero
p_bcr                     = zero
ei_bcr                    = zero
ye_bcr                    = zero
temp_bcr                  = zero
gc_bcr                    = zero
ge_bcr                    = zero
psi0_bcr                  = zero
comp_bcr                  = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_evh1_bound_arrays
                                                                                

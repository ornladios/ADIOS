SUBROUTINE dimension_eos_ls_arrays( nx )
!-----------------------------------------------------------------------
!
!    File:         dimension_eos_ls_arrays
!    Module:       dimension_eos_ls_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the eos ls arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x-array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, eos_ls_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE eos_ls_module
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' EOS LS arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_eos_ls_arrays')

!-----------------------------------------------------------------------
!  Allocate eos_ls_module arrays
!-----------------------------------------------------------------------

ALLOCATE (inpvars(nx,4,2,2,2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'inpvars   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (inpvarsa(nx,4), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'inpvarsa  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize eos_ls_module arrays
!-----------------------------------------------------------------------

inpvars                   = zero
inpvarsa                  = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_eos_ls_arrays
                                                                                

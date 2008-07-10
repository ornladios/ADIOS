SUBROUTINE initialize_global_var
!-----------------------------------------------------------------------
!
!    File:         initialize_global_var
!    Module:       initialize_global_var
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/25/03
!
!    Purpose:
!        To initialize certain global variables before readin.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  edit_module, evh1_global, parallel_module
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : nlog
USE evh1_global, ONLY : ndim, ngeomx, ngeomy, ngeomz, nleftx, nlefty, nleftz, nrightx, nrighty, nrightz
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Global variables have been initialized')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ndim                 = -1
ngeomx               = -1
ngeomy               = -1
ngeomz               = -1
nleftx               = -1
nrightx              = -1
nlefty               = -1
nrighty              = -1
nleftz               = -1
nrightz              = -1

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE initialize_global_var

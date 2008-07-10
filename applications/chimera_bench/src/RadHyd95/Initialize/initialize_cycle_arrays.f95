SUBROUTINE initialize_cycle_arrays
!-----------------------------------------------------------------------
!
!    File:         initialize_cycle_arrays
!    Module:       initialize_cycle_arrays
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
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  cycle_module
!  edit_module, numerical_module, parallel_module
!
!-----------------------------------------------------------------------

USE cycle_module

USE edit_module, ONLY : nlog
USE numerical_module, ONLY : zero
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Cycle variables have been initialized')

!-----------------------------------------------------------------------
!        Initialize cycle_module
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Computational cycles
!-----------------------------------------------------------------------

ncycle                    = 0
ncymax                    = 9000000
nrst                      = 0

!-----------------------------------------------------------------------
!  Hydro subcycling wrt transport
!-----------------------------------------------------------------------

nutrans_trns              = .true.

intnu_trns                = 1
intnur_trns               = 1
ncynu_trns                = 0

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE initialize_cycle_arrays

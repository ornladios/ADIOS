SUBROUTINE dimension_nu_energy_grid_arrays( nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         dimension_nu_energy_grid_arrays
!    Module:       dimension_nu_energy_grid_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the neutrino energy grid arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nez       : neutrino energy array dimension
!  nnu       : neutrino flavor array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, nu_energy_grid_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE nu_energy_grid_module
USE parallel_module, ONLY : myid

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Neutrino energy arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_nu_energy_grid_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!           \\\\\ ALLOCATE NU_ENERGY_GRID_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino group energies
!-----------------------------------------------------------------------

ALLOCATE (nnugp(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nnugp     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (unui(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unui      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dunui(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dunui     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unubi(nez+1), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unubi     '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!         \\\\\ INITIALIZE NU_ENERGY_GRID_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino group energies
!-----------------------------------------------------------------------

nnugp                     = 0
nnugpmx                   = 0

unui                      = zero
dunui                     = zero
unubi                     = zero

WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_nu_energy_grid_arrays

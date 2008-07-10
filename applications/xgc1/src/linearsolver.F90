#if !defined(NO_PETSC)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     
!     Mark Adams 
!     petsc_xgc.F90
!     13 January 2006
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine petsc_init( ierr )
  !-----[--.----+----.----+----.-----------------------------------------]
  !     initialize solve linear system and solver.
  !     
  !     input:
  !     
  !     output:
  !     ierr: error code
  !     
  !     side effects:
  !     - creates objects in global data in 'linear_solver.h'
  !
  !-----[--.----+----.----+----.-----------------------------------------]

  implicit none 
  integer ierr
  !     #include "include/finclude/petscviewer.h"  
  call PetscInitialize('./.petscrc',ierr)
!  call PetscInitialize(0,ierr)
  ! log initialization
  !-----[--.----+----.----+----.-----------------------------------------]
  !call PetscLogEventRegister(init_log, 'Solver init     ',0,ierr)
  !call PetscLogEventRegister(solve_log, 'Linear Solve    ',0,ierr)
  !call PetscLogEventRegister(aux_log, 'Aux 1           ',0,ierr)
  !call PetscLogEventRegister(aux2_log, 'Aux 2           ',0,ierr)
  !-----[--.----+----.----+----.-----------------------------------------]
  
end subroutine petsc_init

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine petsc_end( ierr )
  !-----[--.----+----.----+----.-----------------------------------------]
  !     cleanup common blocks and modules
  !     
  !     input:
  !     output:
  !     ierr: error code
  !     
  !-----[--.----+----.----+----.-----------------------------------------]
  implicit none
  integer ierr 
  call PetscFinalize(ierr)
end subroutine petsc_end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#endif



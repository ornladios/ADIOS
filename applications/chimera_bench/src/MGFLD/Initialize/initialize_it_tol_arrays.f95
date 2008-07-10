SUBROUTINE initialize_it_tol_arrays
!-----------------------------------------------------------------------
!
!    File:         initialize_it_tol_arrays
!    Module:       initialize_it_tol_arrays
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
!  edit_module, it_tol_module, numerical_module, parallel_module
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : nlog
USE it_tol_module
USE numerical_module, ONLY : zero
USE parallel_module, ONLY : myid

IMPLICIT none

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Iteration and tolerence variables have been initialized')

!-----------------------------------------------------------------------
!        Initialize incrmnt_module arrays
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino transport tolerances
!-----------------------------------------------------------------------

iternu                    = 40
itfail                    = 1

tolnut                    = 1.d-02
tolnuye                   = 1.d-02
tolnupsi                  = 1.d-01
tolpsimin                 = 1.d-01
a_prec                    = 1.d-02

!-----------------------------------------------------------------------
!  Hydrodynamic tolerances
!-----------------------------------------------------------------------

it_hy                     = 40
tolhy                     = 1.d-02

WRITE (nlog,101)

RETURN
END SUBROUTINE initialize_it_tol_arrays

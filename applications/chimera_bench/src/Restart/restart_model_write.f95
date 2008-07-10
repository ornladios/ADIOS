SUBROUTINE restart_model_write( ndump, nx, nnu )
!-----------------------------------------------------------------------
!
!    File:         restart_model_write
!    Module:       restart_model_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/01/05
!
!    Purpose:
!      To direct the writing of the model configuration for restarting a simulation.
!
!    Subprograms called:
!  nuclear_write          : writes out the nuclear keys and abundances
!  model_write            : writes out the model configuration
!
!    Input arguments:
!  ndump                  : unit number to write restart dump
!  nx                     : x-array extent
!  nnu                    : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  edit_module, io_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE edit_module, ONLY : nouttmp
USE io_module, ONLY: model_write_hdf5, model_read_hdf5

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ndump         ! unit number to write temporary restart dumps
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!       \\\\\ WRITE MODEL CONFIGURATION TO A RESTART FILE /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Write model configuration and nuclear abundances to a restart file
!-----------------------------------------------------------------------

IF ( nouttmp /= 6 ) THEN
  CALL model_write( ndump, nx, nnu )
ELSE ! nouttmp = 6
  CALL model_write_hdf5()
END IF ! nouttmp /= 6

RETURN
END SUBROUTINE restart_model_write

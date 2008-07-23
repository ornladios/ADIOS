SUBROUTINE restart_keys_write_adios( ndump, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         restart_keys_write
!    Module:       restart_keys_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/01/05
!
!    Purpose:
!      To direct the writing of the simulation keys for restarting a simulation.
!
!    Subprograms called:
!  array_dimensions_write : writes out the array dimensions
!  radhyd_write           : writes out the radhyd keys
!  e_advct_write          : writes out the neutrino energy advection keys
!  eos_write              : writes out the equation of state keys
!  transport_write        : writes out the neutrino transport keys
!  edit_write             : writes out the edit keys
!  hydro_write            : writes out the hydro keys
!
!    Input arguments:
!  ndump                  : unit number to write restart dump
!  nez                    : neutrino energy array extent
!  nnu                    : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  cycle_module, edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE cycle_module, ONLY : nrst
USE edit_module, ONLY : nrrst, nouttmp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

character (len=128), INTENT(in)  :: ndump         ! group handle to write to
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
integer*8 :: io_type, handle
INTEGER                          :: adios_err     ! ADIOS error flag
 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!         \\\\\ WRITE SIMULATION KEYS TO A RESTART FILE /////
!
!-----------------------------------------------------------------------
!call adios_get_group (io_type, 'restart.keys'//char(0))
!call adios_open (handle, io_type, trim(ndump)//char(0))

CALL adios_open (handle,'restart.keys'//char(0),trim(ndump)//char(0),adios_err)

!-----------------------------------------------------------------------
!  Write array_dimensions
!-----------------------------------------------------------------------

CALL array_dimensions_write_adios( handle )

!-----------------------------------------------------------------------
!  Write radyhd keys
!-----------------------------------------------------------------------

CALL radhyd_write_adios( handle )

!-----------------------------------------------------------------------
!  Write neutrino energy advection keys
!-----------------------------------------------------------------------

CALL e_advct_write_adios( handle, nnu )

!-----------------------------------------------------------------------
!  Write equation of state keys
!-----------------------------------------------------------------------

CALL eos_write_adios( handle )

!-----------------------------------------------------------------------
!  Write transport keys
!-----------------------------------------------------------------------

CALL transport_write_adios( handle, nnu )

!-----------------------------------------------------------------------
!  Write edit keys
!-----------------------------------------------------------------------

CALL edit_write_adios( handle, nez, nnu )

!-----------------------------------------------------------------------
!  Write hydro keys
!-----------------------------------------------------------------------

CALL hydro_write_adios( handle )

!-----------------------------------------------------------------------
!  Write muclear keys
!-----------------------------------------------------------------------

CALL nuclear_keys_write_adios( handle )

call adios_close (handle)

RETURN
END SUBROUTINE restart_keys_write_adios

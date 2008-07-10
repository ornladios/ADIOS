SUBROUTINE editng_set( n )
!-----------------------------------------------------------------------
!
!    File:         editng_set
!    Module:       editng_set
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/17/05
!
!    Purpose:
!      To print header for neutrino integral data
!
!    Subprograms call:
!      date_and_time_print
!
!    Input arguments:
!  n          : neutrino flavor
!  nprint     : unit number for print
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!      none
!
!    Include files:
!  edit_module
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : nprint, head

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: n             ! neutrino flavor

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=53), DIMENSION(4) :: ng_header
CHARACTER (len=55), DIMENSION(4) :: ung_header

LOGICAL                          :: first = .true.

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1x)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
  101 FORMAT (13x,a53)
  103 FORMAT (12x,a55/)

!-----------------------------------------------------------------------
!
!                  \\\\\ INITIALIZE TEXT /////
!
!-----------------------------------------------------------------------

IF ( first ) THEN
  ng_header(1)      = 'e-type neutrino data (integral over energy zones)'
  ung_header(1)     = '---------------------------------------------------'
  ng_header(2)      = 'e-type antineutrino data (integral over energy zones)'
  ung_header(2)     = '-------------------------------------------------------'
  ng_header(3)      = 'x-type neutrino data (integral over energy zones)'
  ung_header(3)     = '---------------------------------------------------'
  ng_header(4)      = 'x-type antineutrino data (integral over energy zones)'
  ung_header(4)     = '-------------------------------------------------------'
  first             = .false.
END IF

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( nprint )
WRITE (nprint,101) ng_header(n)
WRITE (nprint,103) ung_header(n)

RETURN
END SUBROUTINE editng_set

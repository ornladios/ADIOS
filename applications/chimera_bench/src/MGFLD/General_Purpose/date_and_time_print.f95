SUBROUTINE date_and_time_print(n_unit)
!-----------------------------------------------------------------------
!
!    File:         date_and_time_print
!    Module:       date_and_time_print
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/05/02
!
!    Purpose:
!      To print the cycle number, date, and time to unit number n_unit.
!
!    Subprograms called:
!      date_and_time
!
!    Input arguments:
!      none
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!      ncycle
!
!    Modules:
!      cycle_module
!
!-----------------------------------------------------------------------

USE cycle_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER(LEN=8)                :: date
CHARACTER(LEN=10)               :: clock,tzone

INTEGER, DIMENSION(9)           :: tvalues
INTEGER                         :: n_unit          ! print unit number

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (/,1x,'cycle =',i8,5x,'date: ',a2,'/',a2,'/',a4,5x,'time: ',a2,':',a2,':',a2/)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Get date and time
!-----------------------------------------------------------------------

CALL date_and_time( date, clock, tzone, tvalues )

!-----------------------------------------------------------------------
!  Print cycle, date, and time
!-----------------------------------------------------------------------

WRITE (n_unit,1001) ncycle,date(5:6),date(7:8),date(1:4),clock(1:2),clock(3:4),clock(5:6)

RETURN
END SUBROUTINE date_and_time_print

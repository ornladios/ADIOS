SUBROUTINE power_SASI_plot( time, t_tb, p_SASI, istat_n, nprint )
!-----------------------------------------------------------------------
!
!    File:         power_SASI_plot
!    Module:       power_SASI_plot
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/08/07
!
!    Purpose:
!      To create files of Legendre modes of the power in convection and
!       the SASI
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  time         : elapsed time
!  t_tb         : time from core bounce
!  p_SASI       : power in the convection or SASI as a function of the Legendre mode
!  istat_n      : open new file flag
!  nprint       : unit numberto print data
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module
!  edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE edit_module, ONLY: head


IMPLICIT NONE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER                               :: nleg = 10     ! highest Legendre polynomial to be used 

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                              :: istat_n       ! open new file flag (passed to subroutine)
INTEGER, INTENT(in)                              :: nprint        ! unit numberto print data

REAL(KIND=double), INTENT(in)                    :: time          ! elapsed time
REAL(KIND=double), INTENT(in)                    :: t_tb          ! time from core bounce
REAL(KIND=double), INTENT(in), DIMENSION(nleg)   :: p_SASI        ! power in the convection or SASI as a function of the Legendre mode

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                                          :: l             ! Legendre polynomial order

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*')/)
  101 FORMAT ('     time         t_bounce   ',10(4x,i3,4x)/)
  103 FORMAT (2es15.8,10(es11.3))

!-----------------------------------------------------------------------
!  Write header if file opened for the first time
!-----------------------------------------------------------------------

IF ( istat_n == 0 ) THEN
  WRITE (nprint,3) head
  WRITE (nprint,5)
  WRITE (nprint,101) (l,l=1,nleg)
END IF ! istat_n == 0
 
!-----------------------------------------------------------------------
!  Write convection and SASI power modes and return
!-----------------------------------------------------------------------
      
WRITE (nprint,103) time, t_tb, (p_SASI(l),l=1,nleg)

RETURN
END SUBROUTINE power_SASI_plot

SUBROUTINE radhyd_read( c_radhyd_data, i_radhyd_data, d_radhyd_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         radhyd_read
!    Module:       radhyd_read
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/02/04
!
!    Purpose:
!      To open and close the files for reading the radhyd keys defining the dimensions,
!       geometry, and the principal parameters of the problem.
!
!
!    Subprograms called:
!  read_pack_radhyd_keys : reads and packs data from file radhyd_keys.d
!
!    Input arguments:
!  nrst                  : cycle number at start or restart
!
!    Output arguments:
!  c_radhyd_data         : character array of radhyd keys
!  i_radhyd_data         : integer array of radhyd keys
!  d_radhyd_data         : 64 bit real array of radhyd keys
!
!    Include files:
!  kind_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module

USE edit_module, ONLY : nprint, nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(out), DIMENSION(20) :: c_radhyd_data ! character array of radhyd keys

INTEGER, INTENT(out), DIMENSION(50)           :: i_radhyd_data ! integer array of radhyd keys

REAL(KIND=double), INTENT(out), DIMENSION(50) :: d_radhyd_data ! 64 bit real array of radhyd keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: nread         ! unit number from which to read transport keys
INTEGER                          :: iskipp        ! echo transport keys read flag
INTEGER                          :: istat         ! open-close file flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' Radhyd keys have been read, initialized, and packed')
 1001 FORMAT (' Error in closing radhyd_keys.d in subroutine radhyd_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nread              = 15
iskipp             = 0

!-----------------------------------------------------------------------
!
!                   \\\\\ READ RADHYD KEYS /////
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/radhyd_keys.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/radhyd_keys.d',STATUS='old')

CALL read_pack_radhyd_keys( nread, nprint, iskipp, c_radhyd_data, i_radhyd_data, &
& d_radhyd_data, nrst )

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)

!-----------------------------------------------------------------------
!  Write diagnostic and stop if problem closing file
!-----------------------------------------------------------------------

IF ( istat /= 0 ) THEN
  WRITE (nlog,1001)
  WRITE (nprint,1001)
  STOP
END IF

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE radhyd_read

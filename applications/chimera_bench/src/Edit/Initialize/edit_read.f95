SUBROUTINE edit_read( i_edit_data, d_edit_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         edit_read
!    Module:       edit_read
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the edit keys prescribing the frequency of the various edits.
!
!
!    Subprograms called:
!  read_pack_edit_keys : reads and packs the edit keys from file edit_keys.d
!
!    Input arguments:
!  nrst                : cycle number at start or restart
!
!    Output arguments:
!  i_edit_data         : integer array of edit keys
!  d_edit_data         : real*8 array of edit keys
!
!    Include files:
!  kind_module, array_module
!  edit_module, edit_keys.d
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nez, nnu

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

INTEGER, INTENT(out), DIMENSION(1200+3*40*nnu) :: i_edit_data  ! integer array of edit keys

REAL(KIND=double), INTENT(out), DIMENSION(50)  :: d_edit_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: nread         ! unit number from which to read transport keys
INTEGER                          :: iskipp        ! echo transport keys read flag
INTEGER                          :: istat         ! open-close file flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' Edit keys have been read, initialized, and packed')
 1001 FORMAT (' Error in closing edit_keys.d in subroutine edit_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nread              = 15
iskipp             = 0

!-----------------------------------------------------------------------
!
!                \\\\\ READ AND PACK EDIT KEYS /////
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/edit_keys.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/edit_keys.d',STATUS='old')

CALL read_pack_edit_keys( nread, nprint, iskipp, nez, nnu, i_edit_data, &
& d_edit_data, nrst )

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001)

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE edit_read

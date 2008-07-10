SUBROUTINE eos_read( c_eos_data, d_eos_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         eos_read
!    Module:       eos_read
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/02/04
!
!    Purpose:
!      To set up the files in order to read in the equation of state keys
!       defining the eos table gridding and initial values of eos variables.
!
!
!    Subprograms called:
!  read_pack_eos_keys : reads and packs eos keys from file eos_keys.d
!
!    Input arguments:
!  nrst               : cycle number at start or restart
!
!    Output arguments:
!  c_eos_data         : character array of eos keys
!  d_eos_data         : 64 bit real array of eos keys
!
!    Include files:
!  kind_module
!  edit_module, parallel_modules
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

CHARACTER (len=1), INTENT(out), DIMENSION(1)  :: c_eos_data  ! character array of eos keys

REAL(KIND=double), INTENT(out), DIMENSION(14) :: d_eos_data  ! 64 bit real array of eos keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: nread         ! unit number from which to read eos keys
INTEGER                          :: iskipp        ! echo eos keys read flag
INTEGER                          :: istat         ! open-close file flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' EOS keys have been read, initialized, and packed')
 1001 FORMAT (' Error in closing hydro_keys.d in subroutine hydro_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nread              = 15
iskipp             = 0

!-----------------------------------------------------------------------
!
!                \\\\\ READ AND PACK EOS KEYS /////
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/eos_keys.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/eos_keys.d',STATUS='old')

CALL read_pack_eos_keys( nread, nprint, iskipp, c_eos_data, d_eos_data, nrst )

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
IF ( istat /= 0 ) WRITE (nlog,1001)

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE eos_read

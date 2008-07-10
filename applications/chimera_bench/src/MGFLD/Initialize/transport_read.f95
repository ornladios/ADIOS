SUBROUTINE transport_read( i_trans_data, d_trans_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         transport_read
!    Module:       transport_read
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the transport keys defining the transport sources, solution
!       tolerances, and time step criteria.
!
!
!    Subprograms called:
!  read_pack_transport_keys : reads and packs the transport keys from
!   file transport_keys.d
!
!    Input arguments:
!  nrst                     : cycle number at start or restart
!
!    Output arguments:
!  i_trans_data             : integer array of transport keys
!  d_trans_data             : real*8 array of transport keys
!
!    Include files:
!  kind_module, array_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nez, nezp1, nnu

USE edit_module, ONLY : nprint, nlog
USE parallel_module, ONLY : myid

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(40+2*nnu)                      :: i_trans_data  ! integer array of transport keys

REAL(KIND=double), INTENT(out), DIMENSION((110+3*nnu+3*nez+1)) :: d_trans_data  ! 64 bit real array of transport keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: nread         ! unit number from which to read transport keys
INTEGER                          :: iskipp        ! echo transport keys read flag
INTEGER                          :: istat         ! open-close file flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' Transport keys have been read, initialized, and packed')
 1001 FORMAT (' Error in closing transport_keys.d in subroutine transport_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nread              = 15
iskipp             = 0

!-----------------------------------------------------------------------
!
!              \\\\\ READ AND PACK TRANSPORT KEYS /////
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/transport_keys.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/transport_keys.d',STATUS='old')

CALL read_pack_transport_keys( nread, nprint, iskipp, nez, nezp1, nnu, i_trans_data, &
& d_trans_data, nrst )

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001)

WRITE (nlog,101)

RETURN
END SUBROUTINE transport_read

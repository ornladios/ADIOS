SUBROUTINE hydro_read( i_hydro_data, d_hydro_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         hydro_read
!    Module:       hydro_read
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
!  read_hydro_keys : reads and packs the hydro keys from hydro_keys.d
!
!    Input arguments:
!  nrst            : cycle number at start or restart
!
!    Output arguments:
!  i_hydro_data    : integer array of hydro keys
!  d_hydro_data    : real*8 array of hydro keys
!
!    Include files:
!  kind_module, array_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx

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

INTEGER, INTENT(out), DIMENSION(30)                :: i_hydro_data  ! integer array of transport keys

REAL(KIND=double), INTENT(out), DIMENSION((30+nx)) :: d_hydro_data  ! 64 bit real array of transport keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: nread         ! unit number from which to read transport keys
INTEGER                          :: iskip         ! echo transport keys read flag
INTEGER                          :: istat         ! open-close file flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' Hydro keys have been read, initialized, and packed')
 1001 FORMAT (' Error in closing hydro_keys.d in subroutine hydro_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nread              = 15
iskip              = 0

!-----------------------------------------------------------------------
!
!                \\\\\ READ AND PACK HYDRO KEYS /////
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/hydro_keys.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/hydro_keys.d',STATUS='old')

CALL read_pack_hydro_keys( nread, nprint, iskip, nx, i_hydro_data, d_hydro_data, nrst )

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001)

IF ( myid == 0 ) WRITE (nlog,101)
RETURN
END SUBROUTINE hydro_read

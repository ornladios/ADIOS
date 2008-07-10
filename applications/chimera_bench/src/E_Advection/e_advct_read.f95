SUBROUTINE e_advct_read( i_e_advct_data, d_e_advct_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         e_advct_read
!    Module:       e_advct_read
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the e_advection keys defining the neutrino energy advection
!       solution, tolerances, and time step criteria.
!
!
!    Subprograms called:
!  read_pack_e_advct_keys : reads and packs the energy advection keys
!   from file e_advct_keys.d
!
!    Input arguments:
!  nrst                   : cycle number at start or restart
!
!    Output arguments:
!  i_e_advct_data         : integer array of e_advct keys
!  d_e_advct_data         : real*8 array of e_advct keys
!
!    Include files:
!  kind_module, array_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nnu

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

INTEGER, INTENT(out), DIMENSION(5)                 :: i_e_advct_data  ! integer array of e_advect keys

REAL(KIND=double), INTENT(out), DIMENSION(5+2*nnu) :: d_e_advct_data  ! 64 bit real array of e_advect keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: nread         ! unit number from which to read transport keys
INTEGER                          :: iskipp        ! echo transport keys read flag
INTEGER                          :: istat         ! open-close file flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' E_advection keys have been read, initialized, and packed')
 1001 FORMAT (' Error in closing e_advct_keys.d in subroutine e_advct_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nread              = 15
iskipp             = 0

!-----------------------------------------------------------------------
!
!         \\\\\ READ AND PACK THE ENERGY ADVECTION KEYS /////
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/e_advct_keys.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/e_advct_keys.d',STATUS='old')

CALL read_pack_e_advct_keys( nread, nprint, iskipp, nnu, i_e_advct_data, &
& d_e_advct_data, nrst )

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001)

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE e_advct_read

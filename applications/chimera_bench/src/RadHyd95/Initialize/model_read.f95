SUBROUTINE model_read( i_model_data, d_model_data1, d_model_data2,   &
& d_model_data3, d_psi_data2, nrst )
!-----------------------------------------------------------------------
!
!    File:         model_read
!    Module:       model_read
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the initial model configuration.
!
!
!    Subprograms called:
!  read_pack_initial_model : reads and packs the model configuration from
!   file initial_model.d
!
!    Input arguments:
!  nrst                    : cycle number at start or restart
!
!    Output arguments:
!  i_model_data            : integer array of initial model data
!  d_model_data1           : 64 bit real array of initial model data
!  d_model_data2           : 64 bit real array of initial model data
!  d_model_data3           : 64 bit real array of initial model data
!  d_psi_data2             : 64 bit real array of initial model data
!
!    Include files:
!  kind_module, array_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nez, nnu, ij_ray_dim, ik_ray_dim

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

INTEGER, INTENT(out), DIMENSION(2)                         :: i_model_data  ! integer array of initial model data

REAL(KIND=double), INTENT(out), DIMENSION(7,nx,ij_ray_dim,ik_ray_dim)         :: d_model_data1 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx)                               :: d_model_data2 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx+1)                             :: d_model_data3 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(out), DIMENSION(2,nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: d_psi_data2   ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: nread         ! unit number from which to read transport keys
INTEGER                          :: iskipp        ! echo transport keys read flag
INTEGER                          :: istat         ! open-close file flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' Model parameters have been read, initialized, and packed')
 1001 FORMAT (' Error in closing initial_model.d in subroutine model_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Initialize.....................................................

nread              = 15
iskipp             = 0

!-----------------------------------------------------------------------
!
!          \\\\\ READ AND PACK MODEL CONFIGURATION DATA /////
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/initial_model.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/initial_model.d',STATUS='old')

CALL read_pack_initial_model( nread, nprint, iskipp, nx, nez, nnu, ij_ray_dim, &
& ik_ray_dim, i_model_data, d_model_data1, d_model_data2, d_model_data3,       &
& d_psi_data2, nrst )

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
  IF ( istat /= 0 ) WRITE (nlog,1001)

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE model_read


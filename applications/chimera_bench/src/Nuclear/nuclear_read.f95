SUBROUTINE nuclear_read( c_nuc_data, i_nuc_data, d_nuc_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         nuclear_read
!    Module:       nuclear_read
!    Type:         Program
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the nuclear keys defining the nuclear sources, solution
!       tolerances, and time step criteria.
!
!
!    Subprograms called:
!  read_pack_initial_nuclear_data : reads and packs nuclear keys and data
!   from file nuclear_keys.d
!
!    Input arguments:
!  nrst                   : cycle number at start or restart
!
!    Output arguments:
!  c_nuc_data             : character array of nuclei
!  i_nuc_data             : integer array of nuclear keys
!  d_nuc_data             : real*8 array of nuclear keys
!
!    Include files:
!  kind_module, array_module
!  edit_module. parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nnc, ij_ray_dim, ik_ray_dim

USE edit_module, ONLY : nprint, nlog
USE parallel_module, ONLY : myid, ierr

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER (len=5), INTENT(out), DIMENSION(nnc)                 :: c_nuc_data  ! character array of nuclei

INTEGER, INTENT(out), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim)                         :: i_nuc_data  ! integer array of edit keys

REAL(KIND=double), INTENT(out), DIMENSION(10+4*nnc+(nnc+4)*nx,ij_ray_dim,ik_ray_dim) :: d_nuc_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: nread         ! unit number from which to read transport keys
INTEGER                          :: iskipp        ! echo transport keys read flag
INTEGER                          :: istat         ! open-close file flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 101 FORMAT (' Nuclear keys have been read, initialized, and packed')
 1001 FORMAT (' Error in closing nuclear_keys.d in subroutine nuc_read')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nread              = 15
iskipp             = 0

!-----------------------------------------------------------------------
!
!              \\\\\ READ AND PACK THE NUCLEAR KEYS /////
!
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/nuclear_keys.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/nuclear_keys.d',STATUS='old')

CALL read_pack_initial_nuclear_data( nread, nprint, iskipp, nx, nnc, &
& ij_ray_dim, ik_ray_dim, c_nuc_data, i_nuc_data, d_nuc_data, nrst )

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
END SUBROUTINE nuclear_read

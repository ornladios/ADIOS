SUBROUTINE read_pack_array_dimensions( c_path_data, n_dim_data )
!-----------------------------------------------------------------------
!
!    File:         read_pack_array_dimensions
!    Module:       read_pack_array_dimensions
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/16/04
!
!    Purpose:
!        To initialize the dimenisons of the arrays used in RadHyd.
!
!    Subprograms called:
!  read_array_dimensons  : reads in array dimension from file
!                           array_dimensions.d
!  read_pack_init        : reads in cycle number to determine how
!                           initialize
!  read_model_dimensions : reads in models dimensions for consistency
!                           checks
!
!    Input arguments:
!        none
!
!    Output arguments:
!  c_path_data           : character array of data path
!  n_dim_data            : array dimension data
!
!    Include files:
!        none
!     
!-----------------------------------------------------------------------

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER(len=128), INTENT(out), DIMENSION(3) :: c_path_data  ! character array of data path

INTEGER, INTENT(out), DIMENSION(20)           :: n_dim_data   ! integer array of array dimensions

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len = 128)               :: rst_keys_file ! character string containing name of restart_keys file
CHARACTER (len = 128)               :: data_path     ! path to the output data directories
CHARACTER (len = 128)               :: log_path      ! path to the simulation log (blank writes to the screen)
CHARACTER (len = 128)               :: reset_path    ! path to write the restart key file, reset.d

CHARACTER (len = 128), DIMENSION(1) :: c_init_data   ! character array of initial data

INTEGER, DIMENSION(2)               :: i_init_data   ! integer array of initial data

INTEGER                             :: nread         ! unit number from which to read transport keys
INTEGER                             :: nprint        ! unit number to print diagnostics
INTEGER                             :: iskipp        ! echo transport keys read flag
INTEGER                             :: istat         ! open-close file flag
INTEGER                             :: nx            ! x-array extent
INTEGER                             :: ny            ! y-array extent
INTEGER                             :: nz            ! z-array extent
INTEGER                             :: nez           ! neutrino energy array extent
INTEGER                             :: nnu           ! neutrino flavor array extent
INTEGER                             :: nnc           ! composition array extent
INTEGER                             :: n_proc        ! number of processors assigned to the run
INTEGER                             :: n_proc_y      ! number of processors assigned to the y-zones
INTEGER                             :: n_proc_z      ! number of processors assigned to the z-zones
INTEGER                             :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER                             :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER                             :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER                             :: k_ray_dim     ! number of radial zones on a processor after swapping with z

INTEGER                             :: nrst          ! cycle number to start simulation
INTEGER                             :: nouttmp       ! dummy variable

INTEGER                             :: imin          ! inner physical x (radial) index
INTEGER                             :: imax          ! outer physical x (radial) index
INTEGER                             :: jmin          ! inner physical y (angular) index
INTEGER                             :: jmax          ! outer physical y (angular) index
INTEGER                             :: kmin          ! inner physical z index
INTEGER                             :: kmax          ! outer physical z index

 1001 FORMAT (' Error in closing array_dimenisons.d in subroutine initialize_array_dimenisons')
 1003 FORMAT (' Error in closing reset.d in subroutine initialize_array_dimenisons')
 1005 FORMAT (' Error in closing Data3/Initial_Data/radhyd_keys.d in subroutine initialize_array_dimenisons')
 1007 FORMAT (' Error in closing Data3/Restart/restart.d in subroutine initialize_array_dimenisons')
 1009 FORMAT (' nouttmp=',i3,' in read_pack_array_dimensions, which is not in the range of 1 - 6; nrst=',i8,' ny=',i3)
 2001 FORMAT (' n_proc =',i6,' < 1')
 2003 FORMAT (' nx,ny,nz,nez,nnu,nnc =',6i6,' one or more of which is < 1')
 2005 FORMAT (' imax-imin+1,n_proc =',2i6,' MOD(imax-imin+1,n_proc) /= 0')
 2007 FORMAT (' jmax-jmin+1,n_proc_y =',2i6,' MOD(jmax-jmin+1,n_proc_y) /= 0')
 2009 FORMAT (' kmax-kmin+1,n_proc_z =',2i6,' MOD(kmax-kmin+1,n_proc_z) /= 0')
 2011 FORMAT (' n_proc_y * n_proc_z,n_proc =',2i6,' MOD(n_proc_y * n_proc_z, n_proc ) /= 0')
 2013 FORMAT (' (jmax-jmin+1) * (kmax-kmin+1),n_proc =',2i6,' n_proc > (jmax-jmin+1) * (kmax-kmin+1)')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Define edit unit number for this read-in opperation
!-----------------------------------------------------------------------

nprint             = 6

!-----------------------------------------------------------------------
!
!                \\\\\ READ IN ARRAY DIMENSIONS /////
!
!-----------------------------------------------------------------------

nread              = 15
iskipp             = 0

!-----------------------------------------------------------------------
!  Open file array_dimensions.d
!-----------------------------------------------------------------------

OPEN (UNIT=nread,FILE='Data3/Initial_Data/array_dimensions.d',STATUS='new',IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/array_dimensions.d',STATUS='old')

!-----------------------------------------------------------------------
!  Read file array_dimensions.d
!-----------------------------------------------------------------------

CALL read_array_dimensions( nread, nprint, iskipp, nx, ny, nz, nez, nnu, &
& nnc, n_proc, n_proc_y, n_proc_z, data_path, log_path, reset_path )

!-----------------------------------------------------------------------
!  Close file array_dimensions.d
!-----------------------------------------------------------------------

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (*,1001)
  WRITE (nprint,1001)
  STOP
END IF ! istat /= 0

!-----------------------------------------------------------------------
!
!                 \\\\\ READ IN MODEL DIMENSIONS /////
!
!  This is a preliminary read to check the consistency of the array
!   dimensions. The model dimensions will read in again later and
!   broadcast to all processors.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Open reset.d
!-----------------------------------------------------------------------

OPEN (UNIT=nread, FILE='Data3/Initial_Data/reset.d', STATUS='new', IOSTAT=istat)
IF ( istat /= 0 ) OPEN (UNIT=nread, FILE='Data3/Initial_Data/reset.d', STATUS='old')

!-----------------------------------------------------------------------
!  Read nrst
!-----------------------------------------------------------------------

CALL read_pack_init( nread, c_init_data, i_init_data, nrst, nouttmp )

!-----------------------------------------------------------------------
!  Close reset.d
!-----------------------------------------------------------------------

CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
IF ( istat /= 0 ) THEN
  WRITE (*,1003)
  WRITE (nprint,1003)
  STOP
END IF ! istat /= 0

IF ( nrst == 0 ) THEN

!-----------------------------------------------------------------------
!  If nrst = 0, open radhyd_keys.d 
!-----------------------------------------------------------------------

  OPEN (UNIT=nread,FILE='Data3/Initial_Data/radhyd_keys.d',STATUS='new',IOSTAT=istat)
  IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/radhyd_keys.d',STATUS='old')

!-----------------------------------------------------------------------
!  Read problem dimensions
!-----------------------------------------------------------------------

  CALL read_model_dimensions( nread, nprint, iskipp, imin, imax, jmin, jmax, &
& kmin, kmax )

!-----------------------------------------------------------------------
!  Close radhyd_keys.d
!-----------------------------------------------------------------------

  CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (*,1005)
    WRITE (nprint,1005)
    STOP
  END IF ! istat /= 0

ELSE ! nrst > 0

  IF ( ny == 1  .and.  nz == 1 ) THEN

!-----------------------------------------------------------------------
!  If nrst > 0 and ny = nz = 1, open radhyd_keys.d
!-----------------------------------------------------------------------

    OPEN (UNIT=nread,FILE='Data3/Initial_Data/radhyd_keys.d',STATUS='new',IOSTAT=istat)
    IF ( istat /= 0 ) OPEN (UNIT=nread,FILE='Data3/Initial_Data/radhyd_keys.d',STATUS='old')
  
  ELSE ! ny > 1 or nz > 1

!-----------------------------------------------------------------------
!  If nrst > 0 and ny /= and/or nz /= 1, open
!    
!     rst_tmp1_keys.d        if nouttmp = 1
!     rst_tmp2_keys.d        if nouttmp = 2
!     restart_keysxxxxxxx.d  if nouttmp = 3
!     restart_final_keys.d'  if nouttmp = 4
!     restart_keysxxxxxxx.d  if nouttmp = 5
!     restart_keysxxxxxxx.d  if nouttmp = 6
!-----------------------------------------------------------------------

    IF ( nouttmp == 1 ) THEN
      OPEN (UNIT=nread, FILE=TRIM(data_path)//'/Restart/rst_tmp1_keys.d', STATUS='new', &
&      IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=nread, FILE=TRIM(data_path)//'/Restart/rst_tmp1_keys.d', &
&      STATUS='old')

    ELSE IF ( nouttmp == 2 ) THEN
      OPEN (UNIT=nread,FILE=TRIM(data_path)//'/Restart/rst_tmp2_keys.d', STATUS='new', &
&      IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=nread, FILE=TRIM(data_path)//'/Restart/rst_tmp2_keys.d', &
&      STATUS='old')
    
    ELSE IF ( nouttmp == 3  .or.  nouttmp == 5  .or.  nouttmp == 6 ) THEN
      WRITE (rst_keys_file,'(a21,i7.7,a2)') '/Restart/restart_keys', &
&      nrst,'.d'
      rst_keys_file     = TRIM(data_path)//TRIM(rst_keys_file)
      OPEN (UNIT=nread, FILE=TRIM(rst_keys_file), STATUS='new', IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=nread, FILE=TRIM(rst_keys_file), STATUS='old')
    
    ELSE IF ( nouttmp == 4 ) THEN
      OPEN (UNIT=nread,FILE=TRIM(data_path)//'/Restart/restart_final_keys.d', STATUS='new', &
&      IOSTAT=istat)
      IF ( istat /= 0 ) OPEN (UNIT=nread, FILE=TRIM(data_path)//'/Restart/restart_final_keys.d', &
&      STATUS='old')

    ELSE
      WRITE (nprint,1009) nouttmp,nrst,ny
      WRITE (*,1009) nouttmp,nrst,ny
      STOP

    END IF ! nouttmp == 1

  END IF ! ny == 1  .and.  nz == 1

!-----------------------------------------------------------------------
!  Read problem dimensions
!-----------------------------------------------------------------------

  CALL read_model_dimensions( nread, nprint, iskipp, imin, imax, jmin, jmax, &
&  kmin, kmax )

!-----------------------------------------------------------------------
!  Close file
!-----------------------------------------------------------------------

  CLOSE (UNIT=nread,STATUS='keep',IOSTAT=istat)
  IF ( istat /= 0 ) THEN
    WRITE (*,1005)
    WRITE (nprint,1005)
    STOP
  END IF ! istat /= 0

END IF ! nrst == 0

!-----------------------------------------------------------------------
!
!           \\\\\ CHECK ARRAY EXTENTS FOR CONSISTENCY /////
!
!-----------------------------------------------------------------------

IF ( n_proc < 1 ) THEN
  WRITE (*,2001) n_proc
  WRITE (nprint,2001) n_proc
  STOP
END IF ! n_proc <= 0

IF ( nx < 1  .or.  ny < 1  .or.  nz < 1  .or.  nez < 1  .or.  nnu < 1  .or.  nnc < 1) THEN
  WRITE (*,2003) nx, ny, nz, nez, nnu, nnc
  WRITE (nprint,2003) nx, ny, nz, nez, nnu, nnc
  STOP
END IF ! nx <= 0, etc.

IF ( MOD( imax-imin+1, n_proc ) /= 0  .and.  MOD( n_proc, imax-imin+1 ) /= 0 ) THEN
  WRITE (*,2005) imax-imin+1, n_proc
  WRITE (nprint,2005) imax-imin+1, n_proc
  STOP
END IF ! MOD( imax-imin+1, n_proc ) /= 0

IF ( MOD( jmax-jmin+1, n_proc_y ) /= 0 ) THEN
  WRITE (*,2007) jmax-jmin+1, n_proc_y
  WRITE (nprint,2007) jmax-jmin+1, n_proc_y
  STOP
END IF ! MOD( jmax-jmin+1, n_proc_y ) /= 0

IF ( MOD( kmax-kmin+1, n_proc_z ) /= 0 ) THEN
  WRITE (*,2009) kmax-kmin+1, n_proc_z
  WRITE (nprint,2009) kmax-kmin+1, n_proc_z
  STOP
END IF ! MOD( kmax-kmin+1, n_proc_z ) /= 0

IF ( MOD( n_proc_y * n_proc_z, n_proc ) /= 0 ) THEN
  WRITE (*,2011) n_proc_y * n_proc_z, n_proc
  WRITE (nprint,2011) n_proc_y * n_proc_z, n_proc
  STOP
END IF ! MOD( n_proc_y * n_proc_z, n_proc ) /= 0

IF ( (jmax-jmin+1) * (kmax-kmin+1) < n_proc ) THEN
  WRITE (*,2013) (jmax-jmin+1) * (kmax-kmin+1), n_proc
  WRITE (nprint,2013) (jmax-jmin+1) * (kmax-kmin+1), n_proc
  STOP
END IF ! (jmax-jmin+1) * (kmax-kmin+1) < n_proc )

!-----------------------------------------------------------------------
!
! \\\\\ COMPUTE I_RAY_DIMY, I_RAY_DIMZ, J_RAY_DIM, AND K_RAY_DIM /////
!
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  j_ray_dim  : number of radial zones on a processor after swapping
!                with y
!  k_ray_dim  : number of radial zones on a processor after swapping
!                with z
!-----------------------------------------------------------------------

ij_ray_dim           = ( jmax - jmin + 1 )/n_proc_y
ik_ray_dim           = ( kmax - kmin + 1 )/n_proc_z
j_ray_dim            = ( imax - imin + 1 )/n_proc_y
k_ray_dim            = ( imax - imin + 1 )/n_proc_z

!-----------------------------------------------------------------------
!
!                  \\\\\ PACK ARRAY DIMENSIONS /////
!
!-----------------------------------------------------------------------

c_path_data(1)     = data_path
c_path_data(2)     = log_path
c_path_data(3)     = reset_path

n_dim_data(1)      = nx
n_dim_data(2)      = ny
n_dim_data(3)      = nz
n_dim_data(4)      = nez
n_dim_data(5)      = nnu
n_dim_data(6)      = nnc
n_dim_data(7)      = n_proc
n_dim_data(8)      = n_proc_y
n_dim_data(9)      = n_proc_z
n_dim_data(10)     = ij_ray_dim
n_dim_data(11)     = ik_ray_dim
n_dim_data(12)     = j_ray_dim
n_dim_data(13)     = k_ray_dim


RETURN
END SUBROUTINE read_pack_array_dimensions

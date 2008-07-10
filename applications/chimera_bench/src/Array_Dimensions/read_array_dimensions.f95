SUBROUTINE read_array_dimensions( nreadp, nprint, iskip, nx, ny, nz, nez, nnu, &
& nnc, n_proc, n_proc_y, n_proc_z, data_path, log_path, reset_path )
!-----------------------------------------------------------------------
!
!    File:         read_array_dimensions
!    Module:       read_array_dimensions
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/01/04
!
!    Purpose:
!      To read in the dimensions of the arrays used in RadHyd.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp     : unit number from which to read
!  nprint     : unit number from which to print
!  iskip      : print key
!
!    Output arguments:
!  nx         : x-array extent
!  ny         : y-array extent
!  nz         : z-array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  nnc        : composition array extent
!  n_proc     : number of processors assigned to the run
!  n_proc_y   : number of processors assigned to the y-zones
!  n_proc_z   : number of processors assigned to the z-zones
!  data_path  : path to the output data directories
!  log_path   : path to the simulation log
!  reset_path : path to write the restart key file, reset.d
!
!    Include files:
!        none
!
!-----------------------------------------------------------------------

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                :: nreadp        ! unit number to read from
INTEGER, INTENT(in)                :: nprint        ! unit number to print to
INTEGER, INTENT(in)                :: iskip         ! print key

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER (len = 128), INTENT(out) :: data_path     ! path to the output data directories
CHARACTER (len = 128), INTENT(out) :: log_path      ! path to the sinulation log
CHARACTER (len = 128), INTENT(out) :: reset_path    ! path to the sinulation log

INTEGER, INTENT(out)               :: nx            ! x-array extent
INTEGER, INTENT(out)               :: ny            ! y-array extent
INTEGER, INTENT(out)               :: nz            ! z-array extent
INTEGER, INTENT(out)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(out)               :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(out)               :: nnc           ! composition array extent
INTEGER, INTENT(out)               :: n_proc        ! number of processors assigned to the run
INTEGER, INTENT(out)               :: n_proc_y      ! number of processors assigned to the y-zones
INTEGER, INTENT(out)               :: n_proc_z      ! number of processors assigned to the z-zones

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                  :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)                 :: name          ! 16 character variable to read in data name
CHARACTER (len=128)                :: line          ! 120 character variable to read in data line

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (a128)
  105 FORMAT (1x,a128)
  107 FORMAT (10x,a118)
  109 FORMAT (a6,4x,a118)
  111 FORMAT (20x,i10)
  113 FORMAT (1x,a6,14x,i10,42x,a16)
  401 FORMAT (' The following data card was not recognized in subroutine read_array_dimensions')
  403 FORMAT (1x,132a)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

data_path                 = ' '
log_path                  = ' '
reset_path                = ' '

!-----------------------------------------------------------------------
!
!                 \\\\\ READ IN ARRAY DIMENSIONS /////
!
!-----------------------------------------------------------------------

READ: DO

!-----------------------------------------------------------------------
!  Read line
!-----------------------------------------------------------------------

  READ (nreadp,101,END=5000) line

  type               = line(1:6)
  name               = line(73:88)

!-----------------------------------------------------------------------
!  Ignore comment cards in the data stream
!-----------------------------------------------------------------------

  IF ( type == 'cccccc'  .or.  type(1:1) == '!'  .or.  type == '      ') THEN
    IF ( iskip == 0 ) WRITE (nprint,105) line
    CYCLE
  END IF ! type = 'cccccc' or type = '!'

!-----------------------------------------------------------------------
!  Array dimensions
!-----------------------------------------------------------------------

!........nx

  IF ( type == 'nx    ' ) THEN
    READ (line ,111) nx
    IF ( iskip == 0 )  WRITE (nprint,113) type,nx,name
    CYCLE
  END IF ! type = 'nz    '

!........ny

  IF ( type == 'ny    ' ) THEN
    READ (line ,111) ny
    IF ( iskip == 0 )  WRITE (nprint,113) type,ny,name
    CYCLE
  END IF ! type = 'ny    '

!........nz

  IF ( type == 'nz    ' ) THEN
    READ (line ,111) nz
    IF ( iskip == 0 )  WRITE (nprint,113) type,nz,name
    CYCLE
  END IF ! type = 'nz    '

!........nez

  IF ( type == 'nez   ' ) THEN
    READ (line ,111) nez
    IF ( iskip == 0 )  WRITE (nprint,113) type,nez,name
    CYCLE
  END IF ! type = 'nez   '

!........nmu

  IF ( type == 'nnu   ' ) THEN
    READ (line ,111) nnu
    IF ( iskip == 0 )  WRITE (nprint,113) type,nnu,name
    CYCLE
  END IF ! type = 'nnu   '

!........nnc

  IF ( type == 'nnc   ' ) THEN
    READ (line ,111) nnc
    IF ( iskip == 0 )  WRITE (nprint,113) type,nnc,name
    CYCLE
  END IF ! type = 'nnc   '

!........n_proc

  IF ( type == 'n_proc' ) THEN

!.............n_proc

    IF ( name == 'n_proc  ') THEN
      READ (line ,111) n_proc
      IF ( iskip == 0 )  WRITE (nprint,113) type,n_proc,name
      CYCLE
    END IF ! name = 'n_proc  '

!.............n_proc_x

    IF ( name == 'n_proc_y') THEN
      READ (line ,111) n_proc_y
      IF ( iskip == 0 )  WRITE (nprint,113) type,n_proc_y,name
      CYCLE
    END IF ! name = 'n_proc_y'

!.............n_proc_y

    IF ( name == 'n_proc_z') THEN
      READ (line ,111) n_proc_z
      IF ( iskip == 0 )  WRITE (nprint,113) type,n_proc_z,name
      CYCLE
    END IF ! name = 'n_proc_z'

  END IF ! type = 'n_proc'

!-----------------------------------------------------------------------
!  data_path
!-----------------------------------------------------------------------

  IF ( type == 'data_p' ) THEN
    READ (line ,107) data_path
    IF ( iskip == 0 ) WRITE (nprint,109) type,data_path
    CYCLE
  END IF ! type = 'data_p'

!-----------------------------------------------------------------------
!  log_path
!-----------------------------------------------------------------------

  IF ( type == 'log_p ' ) THEN
    READ (line ,107) log_path
    IF ( iskip == 0 ) WRITE (nprint,109) type,log_path
    CYCLE
  END IF ! type = 'log_p '

!-----------------------------------------------------------------------
!  reset_path
!-----------------------------------------------------------------------

  IF ( type == 'rset_p' ) THEN
    READ (line ,107) reset_path
    IF ( iskip == 0 ) WRITE (nprint,109) type,reset_path
    CYCLE
  END IF ! type = 'log_p '

!-----------------------------------------------------------------------
!  Unrecognized card
!-----------------------------------------------------------------------

  WRITE (nprint,401)
  WRITE (nprint,403) line
  CYCLE

END DO READ

 5000 RETURN
END SUBROUTINE read_array_dimensions

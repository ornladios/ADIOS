SUBROUTINE read_model_dimensions( nreadp, nprint, iskip, imin, imax, jmin, &
& jmax, kmin, kmax )
!-----------------------------------------------------------------------
!
!    File:         read_model_dimensions
!    Module:       read_model_dimensions
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/02/04
!
!    Purpose:
!      To read in the radhyd keys defining the dimensions, geometry,
!       and principal parameters of the probkem, and pack them in a
!       character array and an integer array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp  : unit number to read from
!  nprint  : unit number to print to
!  iskip   : echo data read flag
!
!    Output arguments:
!  imin    : inner x-array index
!  imax    : outer x-array index
!  jmin    : inner y-array index
!  jmax    : outer y-array index
!  kmin    : inner z-array index
!  kmax    : outer z-array index
!
!    Include files:
!      kind_module
!
!-----------------------------------------------------------------------

USE kind_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out)             :: imin      ! inner x-array index
INTEGER, INTENT(out)             :: imax      ! outer x-array index
INTEGER, INTENT(out)             :: jmin      ! inner y-array index
INTEGER, INTENT(out)             :: jmax      ! outer y-array index
INTEGER, INTENT(out)             :: kmin      ! inner z-array index
INTEGER, INTENT(out)             :: kmax      ! outer z-array index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line

INTEGER                          :: int           ! integer data variable to read in an interger datum
INTEGER                          :: istat         ! read flag

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (a128)
  111 FORMAT (20x,i10)
  113 FORMAT (1x,a6,14x,i10,42x,a16)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                   \\\\\ READ RADHYD DATA /////
!
!-----------------------------------------------------------------------

REWIND (nreadp)

READ: DO

!-----------------------------------------------------------------------
!  Read line
!-----------------------------------------------------------------------

  READ (nreadp,101,IOSTAT=istat) line
  IF ( istat < 0 ) EXIT

  type               = line(1:6)
  name               = line(73:88)

!-----------------------------------------------------------------------
!  Ignore comment cards in the data stream
!-----------------------------------------------------------------------

  IF ( type == 'cccccc'  .or.  type(1:1) == '!'  .or.  type == '      ') CYCLE

!-----------------------------------------------------------------------
!  Spatial dimensions of the problem
!-----------------------------------------------------------------------

!........evh1zn

  IF ( type == 'evh1zn' ) THEN
    READ (line ,111) int

    IF ( name == 'imin    ' ) THEN
      imin           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,imin,name
      CYCLE
    END IF ! name = 'imin    '

    IF ( name == 'imax    ' ) THEN
      imax           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,imax,name
      CYCLE
    END IF ! name = 'imax    '

    IF ( name == 'jmin    ' ) THEN
      jmin           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,jmin,name
      CYCLE
    END IF ! name = 'jmin    '

    IF ( name == 'jmax    ' ) THEN
      jmax           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,jmax,name
      CYCLE
    END IF ! name = 'jmax    '

    IF ( name == 'kmin    ' ) THEN
      kmin           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,kmin,name
      CYCLE
    END IF ! name = 'kmin    '

    IF ( name == 'kmax    ' ) THEN
      kmax           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,kmax,name
      CYCLE
    END IF ! name = 'kmax    '

  END IF ! type = 'evh1zn'

END DO READ

RETURN
END SUBROUTINE read_model_dimensions

SUBROUTINE read_pack_eos_keys( nreadp, nprint, iskip, c_eos_data, d_eos_data, &
& nrst )
!-----------------------------------------------------------------------
!
!    File:         read_pack_eos_keys
!    Module:       read_pack_eos_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the equation of state keys defining the table gridding
!       and to pack these into a character array and a real*8 arrat.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp     : unit number from which to read
!  nprint     : unit number from which to print
!  iskip      : echo data read flag
!  nrst       : cycle number at start or restart
!
!    Output arguments:
!  c_eos_data : character array of eos keys
!  d_eso_data : real*8 array of eos keys
!
!    Include files:
!  kind_module, numerical_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER (len=1), INTENT(out), DIMENSION(1)  :: c_eos_data  ! character array of eos keys

REAL(KIND=double), INTENT(out), DIMENSION(14) :: d_eos_data  ! 64 bit real array of eos keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=1)                :: char1         ! variable for reading in 1 character data
CHARACTER (len=6)                :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line

INTEGER                          :: i             ! do index

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum

!-----------------------------------------------------------------------
!  Equation of state identifier
!-----------------------------------------------------------------------
!  eos_i : equation of state identifier
!
!     eos_i = 'L'  : Lattimer-Swesty equation of state is used.
!     eos_i = 'B'  : Cooperstein-BCK equation of state is used.
!-----------------------------------------------------------------------

CHARACTER (LEN=1)                                       :: eos_i

!-----------------------------------------------------------------------
!  Cube grid
!-----------------------------------------------------------------------
!  dgrid(i), tgrid(i), ygrid(i) : log(rho), log(t), ye space is overlain with a uniform grid of 
!   'dgrid(i)' divisions per unit change in log(rho), 'tgrid(i)' divisions per unit change in 
!   log(t), and 'ygrid(i)' divisions per 0.5 change in ye. Equation of state, nuclear reaction
!   rate, and neutrino interaction variables at each radial zone are interpolated from values at
!   the nearest corners on this grid.
!
!  rhoes(k) : The variables dgrid, tgrid, and ygrid are each 3 element arrays permitting different 
!   partitionings of log(rho), log(t), ye space in different density regimes delimited by rhoes. 
!   These different regimes are
!
!     regime 1:             rho < rhoes(1)
!     regime 2:        rhoes(1) < rho < rhoes(2)
!     regime 3:             rhoes(2) < rho
!
!  idty(j)  : the density regime (i.e., 1, 2, or 3) of radial zone j as given by the above 
!   inequalities.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(3)                         :: dgrid
REAL(KIND=double), DIMENSION(3)                         :: tgrid
REAL(KIND=double), DIMENSION(3)                         :: ygrid
REAL(KIND=double), DIMENSION(2)                         :: rhoes

!-----------------------------------------------------------------------
!  Equation of state borders
!-----------------------------------------------------------------------
!  eosrho : the border density between the LS EOS and the BCK EOS (/fm3).
!-----------------------------------------------------------------------

REAL(KIND=double)                                       :: eosrho

!-----------------------------------------------------------------------
!  NSE flashing and deflashing controls
!-----------------------------------------------------------------------
!  tnse  : temperature at which material is flashed to nse (K)
!
!  tdnse : temperature at which material is deflashed from nse (K)
!-----------------------------------------------------------------------

REAL(KIND=double)                                       :: tnse
REAL(KIND=double)                                       :: tdnse

  101 FORMAT (a128)
  103 FORMAT (a6)
  105 FORMAT (1x,a128)
  131 FORMAT (20x,i10,5x,e15.8)
  133 FORMAT (1x,a6,14x,i10,5x,es15.8,22x,a16)
  141 FORMAT (35x,e15.8)
  143 FORMAT (1x,a6,29x,es15.8,22x,a16)
  191 FORMAT (29x,a1)
  193 FORMAT (1x,a6,14x,9x,a1,42x,a16)
  401 FORMAT (' The following data card was not recognized')
  403 FORMAT (1x,132a)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

dgrid                = zero
tgrid                = zero
ygrid                = zero
rhoes                = zero
eosrho               = zero
tnse                 = zero
tdnse                = zero

!-----------------------------------------------------------------------
!
!                 \\\\\ READ EOS KEYS /////
!
!-----------------------------------------------------------------------

REWIND (nreadp)

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
!
!               \\\\\ EQUATION OF STATE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  esgrid
!-----------------------------------------------------------------------

  IF ( type == 'esgrid' ) THEN

    IF ( name == 'dgrid   ') THEN
       READ (line ,131) i,dgrid(i)
      IF ( iskip == 0 ) WRITE (nprint,133) type,i,dgrid(i),name
      CYCLE
    END IF ! name = 'dgrid   '

    IF ( name == 'tgrid   ') THEN
      READ (line ,131) i,tgrid(i)
      IF ( iskip == 0 ) WRITE (nprint,133) type,i,tgrid(i),name
      CYCLE
    END IF ! name = 'tgrid   '

    IF ( name == 'ygrid   ') THEN
      READ (line ,131) i,ygrid(i)
      IF ( iskip == 0 ) WRITE (nprint,133) type,i,ygrid(i),name
      CYCLE
    END IF ! name = 'ygrid   '

    IF ( name == 'rhoes   ') THEN
      READ (line ,131) i,rhoes(i)
      IF ( iskip == 0 ) WRITE (nprint,133) type,i,rhoes(i),name
      CYCLE
    END IF ! name = 'rhoes   '

  END IF ! type = 'esgrid'

!-----------------------------------------------------------------------
!  eos_i
!-----------------------------------------------------------------------

  IF ( type == 'eos_i ' ) THEN
    READ (line ,191) char1
    eos_i             = char1
    IF ( iskip == 0 ) WRITE (nprint,193) type,eos_i,name
    CYCLE
  END IF ! type = 'eos_i '

!-----------------------------------------------------------------------
!  eosrho
!-----------------------------------------------------------------------

  IF ( type == 'eosrho' ) THEN
    READ (line ,141) rl
    eosrho            = rl
    IF ( iskip == 0 ) WRITE (nprint,143) type,eosrho,name
    CYCLE
  END IF ! type = 'eosrho'

!-----------------------------------------------------------------------
!  tnse
!-----------------------------------------------------------------------

  IF ( type == 'tnse  ' ) THEN
    READ (line ,141) rl
    tnse              = rl
    IF ( iskip == 0 ) WRITE (nprint,143) type,tnse,name
    CYCLE
  END IF ! type = 'tnse  '

!-----------------------------------------------------------------------
!  tdnse
!-----------------------------------------------------------------------

  IF ( type == 'tdnse ' ) THEN
    READ (line ,141) rl
    tdnse             = rl
    IF ( iskip == 0 ) WRITE (nprint,143) type,tdnse,name
    CYCLE
  END IF ! type = 'tdnse '

!-----------------------------------------------------------------------
!  Unrecognized card
!-----------------------------------------------------------------------

  IF ( nrst == 0 ) THEN
    WRITE (nprint,401)
    WRITE (nprint,403) line
  END IF

END DO READ

!-----------------------------------------------------------------------
!
!                 \\\\\ PACK EOS KEYS /////
!
!-----------------------------------------------------------------------

 5000 CONTINUE

c_eos_data(1)                     = eos_i

DO i = 1,3
  d_eos_data(0*3+i)               = dgrid(i)
  d_eos_data(1*3+i)               = tgrid(i)
  d_eos_data(2*3+i)               = ygrid(i)
END DO

d_eos_data(10)                    = rhoes(1)
d_eos_data(11)                    = rhoes(2)
d_eos_data(12)                    = eosrho
d_eos_data(13)                    = tnse
d_eos_data(14)                    = tdnse

RETURN
END SUBROUTINE read_pack_eos_keys

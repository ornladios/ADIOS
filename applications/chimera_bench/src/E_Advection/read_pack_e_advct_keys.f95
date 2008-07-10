SUBROUTINE read_pack_e_advct_keys( nreadp, nprint, iskip, nnu, i_e_advct_data, &
& d_e_advct_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         read_pack_e_advct_keys
!    Module:       read_pack_e_advct_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/26/04
!
!    Purpose:
!      To read in the e_advection keys turning on the neutrino e_advection and
!       determining the tolerances, and to pack them into an integer and a real*8 array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp         : unit number from which to read
!  nprint         : unit number from which to print
!  iskip          : echo data read flag
!  nnu            : neutrino flavor dimension
!  nrst           : cycle number at start or restart
!
!    Output arguments:
!  i_e_advct_data : integer array of e_advct keys
!  d_e_advct_data : real*8 array of e_advct keys
!
!    Include files:
!  kind_module, numerical_module
!  edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor dimension
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(5)                 :: i_e_advct_data  ! integer array of e_advect keys

REAL(KIND=double), INTENT(out), DIMENSION(5+2*nnu) :: d_e_advct_data  ! 64 bit real array of e_advect keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 128 character variable to read in data line
CHARACTER (len=15)               :: var_name

INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum

!-----------------------------------------------------------------------
!  Neutrino energy advection controls
!-----------------------------------------------------------------------
!  ivc_x     : neutrino energy advection due to x-hydro switch.
!
!   ivc_x = 0 : x-hydro neutrino energy advection bypassed.
!   ivc_x = 1 : x-hydro neutrino energy advection included.
!
!  ivc_y     : neutrino energy advection due to y-hydro switch.
!
!   ivc_y = 0 : y-hydro neutrino energy advection bypassed.
!   ivc_y = 1 : y-hydro neutrino energy advection included.
!
!  ivc_z     : neutrino energy advection due to y-hydro switch.
!
!   ivc_z = 0 : z-hydro neutrino energy advection bypassed.
!   ivc_z = 1 : z-hydro neutrino energy advection included.
!
!  tau_advct : optical depth above which neutrinos are advected in
!   energy like a gamma = 4/3 gas
!
!  rhomin_y_eadvect : density below which y-neutrino energy advection
!   is turned off
!
!  rhomin_z_eadvect : density below which z-neutrino energy advection
!   is turned off
!-----------------------------------------------------------------------

INTEGER                                                    :: ivc_x
INTEGER                                                    :: ivc_y
INTEGER                                                    :: ivc_z

REAL(KIND=double)                                          :: tau_advct
REAL(KIND=double)                                          :: rhomin_y_eadvect
REAL(KIND=double)                                          :: rhomin_z_eadvect

!-----------------------------------------------------------------------
!        Time step control criteria
!-----------------------------------------------------------------------
!
!  t_cntrl_e_advct(n)     : n-neutrino zero-moment change time step
!   criterion due to advection, i.e., maximum permitted
!   abs( dpsi0(j,k,n)/psi0(j,k,n) ) due to advection. (n=1,2,3).
!
!  dpsivmx(n), psivmin(n) : parameters used in determining the psi0
!   change time step due to compression, expansion and advection.
!
!     dpsivmx(n) = max( abs( dpsi0(j,k,n) )/( psi0(j,k,n) + psivmin(n) ) ) .
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: t_cntrl_e_advct
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)               :: psivmin

  101 FORMAT (a128)
  105 FORMAT (1x,a128)
  111 FORMAT (20x,i10)
  113 FORMAT (1x,a6,14x,i10,42x,a16)
  131 FORMAT (20x,i10,5x,e15.8)
  133 FORMAT (1x,a6,14x,i10,5x,es15.8,22x,a16)
  141 FORMAT (35x,e15.8)
  143 FORMAT (1x,a6,29x,es15.8,22x,a16)
  401 FORMAT (' The following data card was not recognized')
  403 FORMAT (1x,132a)
 1001 FORMAT (' Allocation problem for array ',a10,' in read_pack_e_advct_keys')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_pack_e_advct_keys')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (t_cntrl_e_advct(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_cntrl_e_advct'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psivmin(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psivmin        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

ivc_x                             = 1
ivc_y                             = 1
ivc_z                             = 1

tau_advct                         = zero
rhomin_y_eadvect                  = 1.d+12
rhomin_z_eadvect                  = 1.d+12
t_cntrl_e_advct                   = zero
psivmin                           = 1.0d-01

i_e_advct_data                    = 0
d_e_advct_data                    = zero

!-----------------------------------------------------------------------
!
!                 \\\\\ READ E_ADVECTION KEYS /////
!
!-----------------------------------------------------------------------

REWIND (nreadp)

READ: DO

!-----------------------------------------------------------------------
!  Read line
!-----------------------------------------------------------------------

  READ (nreadp,101,END=5000) line

  type                            = line(1:6)
  name                            = line(73:88)

!-----------------------------------------------------------------------
!  Ignore comment cards in the data stream
!-----------------------------------------------------------------------

  IF ( type == 'cccccc'  .or.  type(1:1) == '!'  .or.  type == '      ') THEN
    IF ( iskip == 0 ) WRITE (nprint,105) line
    CYCLE
  END IF ! type = 'cccccc' or type = '!'

!-----------------------------------------------------------------------
!  Problem control parameters
!-----------------------------------------------------------------------

!........ivc

  IF ( type == 'ivc   ' ) THEN

!........ivc_x

    IF ( name == 'ivc_x') THEN
      READ (line ,111) ivc_x
      IF ( iskip == 0 )  WRITE (nprint,113) type,ivc_x,name
      CYCLE
    END IF ! type = 'ivc_x   '

!........ivc_y

    IF ( name == 'ivc_y') THEN
      READ (line ,111) ivc_y
      IF ( iskip == 0 )  WRITE (nprint,113) type,ivc_y,name
      CYCLE
    END IF ! type = 'ivc_y   '

!........ivc_z

    IF ( name == 'ivc_z') THEN
      READ (line ,111) ivc_z
      IF ( iskip == 0 )  WRITE (nprint,113) type,ivc_z,name
      CYCLE
    END IF ! type = 'ivc_z   '

  END IF ! type = 'ivc'

!........tau_advct

  IF ( type == 'tau   ' ) THEN
    READ (line ,141) tau_advct
    IF ( iskip == 0 ) WRITE (nprint,143) type,tau_advct,name
    CYCLE
  END IF ! type = 'tau   '

!........rhomin

  IF ( type == 'rhomin' ) THEN

!........rhomin_y_eadvect

    IF ( name == 'rhomin_y' ) THEN
      READ (line ,141) rhomin_y_eadvect
      IF ( iskip == 0 ) WRITE (nprint,143) type,rhomin_y_eadvect,name
      CYCLE
    END IF ! name = 'rhomin_y'

!........rhomin_z_eadvect

    IF ( name == 'rhomin_z' ) THEN
      READ (line ,141) rhomin_z_eadvect
      IF ( iskip == 0 ) WRITE (nprint,143) type,rhomin_z_eadvect,name
      CYCLE
    END IF ! type = 'rhomin_z'

  END IF ! type = 'rhomin'

!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------

!........tcntrl

  IF ( type == 'tcntrl' ) THEN

    READ (line ,131) n,rl

    IF ( name == 't_cntrl_e_advct') THEN
      t_cntrl_e_advct(n) = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,n,t_cntrl_e_advct(n),name
      CYCLE
    END IF ! name = 't_cntrl_e_advct'

  END IF ! type = 'tcntrl'

!........psivtl

  IF ( type == 'psivtl' ) THEN
    READ (line ,131) n,rl
    psivmin(n)  = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,n,psivmin(n),name
    CYCLE
  END IF ! type = 'psivtl'


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
!                 \\\\\ PACK E_ADVECTION KEYS /////
!
!-----------------------------------------------------------------------

 5000 CONTINUE

i_e_advct_data(1)                 = ivc_x
i_e_advct_data(2)                 = ivc_y
i_e_advct_data(3)                 = ivc_z

d_e_advct_data(1)                 = tau_advct
d_e_advct_data(2)                 = rhomin_y_eadvect
d_e_advct_data(3)                 = rhomin_z_eadvect

DO n = 1,nnu
  d_e_advct_data(5+0*nnu+n)       = t_cntrl_e_advct(n)
  d_e_advct_data(5+1*nnu+n)       = psivmin(n)
END DO

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (t_cntrl_e_advct, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 't_cntrl_e_advct'; WRITE (nlog,1001) var_name; END IF
DEALLOCATE (psivmin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psivmin        '; WRITE (nlog,1001) var_name; END IF

RETURN
END SUBROUTINE read_pack_e_advct_keys

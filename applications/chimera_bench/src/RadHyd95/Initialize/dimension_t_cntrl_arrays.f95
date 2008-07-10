SUBROUTINE dimension_t_cntrl_arrays( nx, nnu )
!-----------------------------------------------------------------------
!
!    File:         dimension_t_cntrl_arrays
!    Module:       dimension_t_cntrl_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the time and time step control arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x (radial) array dimension
!  nez       : neutrino energy array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, parallel_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE t_cntrl_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=12)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Radhyd time and time step control arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a12,' in dimension_t_cntrl_arrays')

!-----------------------------------------------------------------------
!        Allocate t_cntrl_module arrays
!-----------------------------------------------------------------------

ALLOCATE (dtj(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtj         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dtjr(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtjr        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dtau_nuadvct(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtau_nuadvct'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dtau_nutrns(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtau_nutrns '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dtauj_nutrns(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtauj_nutrns'; WRITE (nlog,1001) var_name; END IF

ALLOCATE (cdt(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cdt         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cdtinv(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cdtinv      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psimin(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psimin      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psipmin(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psipmin     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsisp(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsisp      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dth(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dth         '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ INITIALIZE T_CNTRL_MODULE_ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Times
!-----------------------------------------------------------------------

time                      = zero
t_start                   = zero
t_bounce                  = zero
t_stop                    = 1.0d+100
tb_stop                   = 1.0d+100

!-----------------------------------------------------------------------
!  Time steps
!-----------------------------------------------------------------------

dtnph                     = 1.d+20
dtnmh                     = 1.d+20
dt_sub                    = 1.d+20
dtnph_trans               = 1.d+20
dtnmhn_trans              = 1.d+20
delt_trans                = 1.d+20
dtj                       = 1.d+20
dtjr                      = 1.d+20
dtau_nuadvct              = 1.d+20
dtau_nutrns               = 1.d+20
dtauj_nutrns              = 1.d+20
dtime_hydro               = 1.d+20
dtime_nuadvct             = 1.d+20
dtime_nutrns              = 1.d+20

!-----------------------------------------------------------------------
!  c x dt
!-----------------------------------------------------------------------

cdt                       = 1.d+20
cdtinv                    = zero

!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------

jrdt                      = 0
jadt                      = 0
tcntrl                    = 1.d+20
dt                        = 1.d+20

!-----------------------------------------------------------------------
!  Time step tolarances
!-----------------------------------------------------------------------

rdtmax                    = 1.2d+00
dtst1                     = zero
dtst2                     = 1.d+20
dtst3                     = 1.d+20
ttst1                     = zero

psimin                    = 1.d-01
psipmin                   = 1.d-01
dpsisp                    = zero

!-----------------------------------------------------------------------
!  Individual time steps
!-----------------------------------------------------------------------

idtj                      = 0
jtv                       = 0
jtvshk                    = 0
ncyrho                    = 0
ncyshk                    = 0

rhojtv                    = 1.d+20
dth                       = 1.d+20

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_t_cntrl_arrays

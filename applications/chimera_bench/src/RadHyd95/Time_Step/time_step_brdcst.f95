SUBROUTINE time_step_brdcst( dt_process, j_radial_dt, j_angular_dt, &
& j_azimuthal_dt, dtnmh, dtnph, dtnph_trans )
!-----------------------------------------------------------------------
!
!    File:         time_step_brdcst
!    Module:       time_step_brdcst
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/20/05
!
!    Purpose:
!      To distribute time step information to all nodes.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  dt_process     : minimum time step for process i
!  j_radial_dt    : radial ray responsible for each dt_process
!  j_angular_dt   : angular ray responsible for each dt_process
!  j_azimuthal_dt : azimuthal ray responsible for each dt_process
!  dtmph          : time step of the preceding cycle
!  dtnph          : time step for the upcoming cycle
!  dtnph_trans    : new source and transport time step
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE t_cntrl_module, ONLY: dtnmh_t=>dtnmh, dtnph_t=>dtnph, dt, jrdt, jadt, &
& jzdt, dtnph_trans_t=>dtnph_trans, dtime_hydro, dtime_nuadvct, dtime_nutrns

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in), DIMENSION(50)           :: j_radial_dt    ! radial ray responsible for each dt_process
INTEGER, INTENT(in), DIMENSION(50)           :: j_angular_dt   ! angular ray responsible for each dt_process
INTEGER, INTENT(in), DIMENSION(50)           :: j_azimuthal_dt ! azimuthal ray responsible for each dt_process

REAL(KIND=double), INTENT(in), DIMENSION(50) :: dt_process     ! time step restricted by all criteria
REAL(KIND=double), INTENT(in)                :: dtnmh          ! time step of the preceding cycle
REAL(KIND=double), INTENT(in)                :: dtnph          ! time step for the upcoming cycle
REAL(KIND=double), INTENT(in)                :: dtnph_trans    ! new source and transport time step

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                      :: i              ! do index

REAL(KIND=double), PARAMETER                 :: dtmax = 1.d+20 ! times step without criteria

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer time step variables to t_cntrl_module
!-----------------------------------------------------------------------

dtnmh_t            = dtnmh
dtnph_t            = dtnph
dtnph_trans_t      = dtnph_trans
dt                 = dt_process
jrdt               = j_radial_dt
jadt               = j_angular_dt
jzdt               = j_azimuthal_dt

!-----------------------------------------------------------------------
!  Deterimine minimum hydro time step
!-----------------------------------------------------------------------

dtime_hydro        = 1.d+20

DO i = 1,10
  IF ( dtime_hydro > dt_process(i) ) dtime_hydro = dt_process(i)
END DO

!-----------------------------------------------------------------------
!  Deterimine minimum source and transport time step
!-----------------------------------------------------------------------

dtime_nutrns       = 1.d+20

DO i = 11,24
  IF ( dtime_nutrns > dt_process(i) ) dtime_nutrns = dt_process(i)
END DO

!-----------------------------------------------------------------------
!  Deterimine minimum neutrino energy advection time step
!-----------------------------------------------------------------------

dtime_nuadvct      = 1.d+20

DO i = 41,48
  IF ( dtime_nuadvct > dt_process(i) ) dtime_nuadvct = dt_process(i)
END DO

RETURN
END SUBROUTINE time_step_brdcst

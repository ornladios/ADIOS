SUBROUTINE hydro_time_step_initialize( jr_min, jr_max, dtime_hydro )
!-----------------------------------------------------------------------
!
!    File:         hydro_time_step_initialize
!    Module:       hydro_time_step_initialize
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To initialize the hydro time step selection.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min      : minimum radial zone index
!  jr_max      : maximum radial zone index
!
!    Output arguments:
!  dtime_hydro : initialized minimum times step by hydro criteria
!
!    Include files:
!  kind_module, numerical_module
!  t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero, one

USE t_cntrl_module, ONLY: dtnmh, dtnph, dtj, dtjr, tcntrl, rdtmax
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input-output variables.
!-----------------------------------------------------------------------

REAL(KIND=double)        :: dtime_hydro     ! minimum times step by hydro criteria

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

LOGICAL                  :: first = .true.

INTEGER                  :: j             ! radial zone index
INTEGER                  :: jr_min          ! minimum zone index
INTEGER                  :: jr_max          ! maximum zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialization at problem start or restart
!-----------------------------------------------------------------------

IF ( first ) THEN

  first            = .false.

  IF ( rdtmax < zero ) rdtmax = one + tcntrl(10)

!........Initialize dtj array

  DO j = jr_min,jr_max
    dtj(j)         = dtnmh
  END DO
  
END IF ! first

!-----------------------------------------------------------------------
!  Initialization during problem execution
!-----------------------------------------------------------------------

!........Switch time step and time step arrays

DO j = jr_min,jr_max
  dtjr(j)          = dtj(j)
END DO

!........Initialize time step

dtime_hydro        = tcntrl(50)

RETURN
END SUBROUTINE hydro_time_step_initialize

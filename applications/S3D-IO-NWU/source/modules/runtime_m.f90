!=========================================================================================
  module runtime_m
!=========================================================================================
! module for run-time variables

  implicit none
!-----------------------------------------------------------------------------------------
! integers

  integer i_time        !time step counter
  integer i_restart     !new run or restart switch: 0 for new run, 1 for restart
  integer i_time_end    !ending time step

  integer i_time_save   !timestep frequency by which to write savefiles
  integer i_write       !switch for writing to screen or file
  integer i_time_mon    !timestep frequency by which to write controller and monitor info
  integer i_time_res    !timestep frequency by which to check resolution
  integer i_time_tec    !timestep frequency by which to write current tecplot file

! reals

  real time             !current time (non-dimensional)
  real time_accum       !current time including substage time (non-dimensional)
  real tstep            !timestep (non-dimensional)
  real time_save        !time at which to write savefiles (seconds)
  real time_save_inc    !increment by which to write savefiles (seconds)
  real time_restart     !time from which the restart was done (seconds)

! characters

  character run_title*20  !unique title of run
!-----------------------------------------------------------------------------------------
  end module runtime_m

!=========================================================================================
  module rk_m
!=========================================================================================
! Module for Explicit Runge-Kutta Time Integration
!
! BUG FIX
! Evatt Hawkes 11-AUG-04
! Divide by zero could occur in timestep calculation in rare cases.
! added parameter vsmall to avoid this.

  implicit none
!-----------------------------------------------------------------------------------------
! set integration type

  character(len=3), parameter :: rk_type ='erk'
!-----------------------------------------------------------------------------------------
! input integers

  integer rk_method       !ERK method (input)
  integer cont_switch     !switch for controller (input)
  integer cfl_switch      !switch for cfl check (input)
  integer i_time_cont     !timestep freqency to write controller info to ts.dat (input)

! other integers

  integer nstage          !number of stages
  integer cont_n_reg      !number of registers in q_err array based on type

! input reals

  real rk_rtol            !RK relative tolerance (input)
  real rk_atol            !RK absolute tolerance (input)
  real k_I                !integral gain (input)
  real k_P                !proportional gain (input)
  real k_D                !derivative gain (input)
  real k_D2               !second derivative gain (input)
  real cont_safety        !factor of safety for controller (input)
  real tstep_init         !initial time step in seconds (input)
  real tstep_min          !minimum time step in seconds (input)
  real tstep_max          !maximum time step in seconds (input)
  real cont_factor        !overall maximum increase in timestep (input)

! other reals

  real rk_p               !order of embedded method
  real cfl_no             !invicid cfl limit
  real fo_no              !viscous cfl limit
  real q_err_max          !maximum Runge-Kutta error
  real rk_tol             !overall RK tolerance (based on rk_atol and rk_rtol)

! real arrays

  real, allocatable :: rk_alpha(:)          !rk integeration coefficients
  real, allocatable :: rk_beta(:)           !rk integeration coefficients
  real, allocatable :: rk_gamma(:)          !rk integeration coefficients
  real, allocatable :: rk_delta(:)          !rk integeration coefficients
  real, allocatable :: rk_eps(:)            !rk integeration coefficients
  real, allocatable :: rk_time(:)           !rk sub-stage time
  real, allocatable :: rk_err(:)            !rk error coefficients
  real, allocatable :: q_err(:,:,:,:,:)     !solution vector error
  real, allocatable :: tstep_vec(:)         !vector of timesteps (current and previous)

! characters

  character*4 cont_type     !controller type: I, PI, PID, PID2, etc...
!-----------------------------------------------------------------------------------------
  end module rk_m

SUBROUTINE parameter_reset( ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         parameter_reset
!    Module:       parameter_reset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/20/05
!
!    Purpose:
!      To reset parameters based on certain criteria during a simulation
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  array_module, numerical_module
!  cycle_module, evh1_global, edit_module, parallel_module, 
!  prb_cntl_module, radial_ray_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE array_module, ONLY : nnu
USE numerical_module, ONLY : zero

USE cycle_module, ONLY : ncycle, intnu_trns
USE evh1_global, ONLY : lagrangian
USE edit_module, ONLY : nprint, intedc, nedc, intede, nede, intdma, nedma, &
& intedh, nedh, intdps, nedps, intedu, nedu, intedy, nedy, intdsc, nedsc, &
& intedn, nedn, intdng, nedng, intprt, nprt, nlog
USE parallel_module, ONLY : myid
USE prb_cntl_module, ONLY : lagr_p=>lagr
USE radial_ray_module, ONLY : t_bounce, lagr, t_bounce_lagr_chg, &
& m_grid, t_bounce_mgrd_chg, time, nedc_r=>nedc, nede_r=>nede, nedmi_r=>nedmi, &
& nedma_r=>nedma, nedh_r=>nedh, nedps_r=>nedps, nedu_r=>nedu, nedy_r=>nedy, &
& nedsc_r=>nedsc, nedn_r=>nedn, nedng_r=>nedng, nu_equil, t_equilibrate
USE t_cntrl_module, ONLY : tcntrl

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray_dim      ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i               ! do index
INTEGER                          :: n               ! neutrino flavor index

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' lagr has been set to "no" in parameter_reset')
  103 FORMAT (' m_grid has been set to "no" in parameter_reset')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Check criteria for switching lagr to 'no'
!-----------------------------------------------------------------------

IF ( lagr == 'ye' ) THEN
  IF ( t_bounce > zero  .and.  time - t_bounce > t_bounce_lagr_chg ) THEN
    lagr                    = 'no'
    lagr_p                  = 'no'
    lagrangian              = .false.
    IF ( myid == 0 ) WRITE (nprint,101)
    WRITE (nlog,101)
  END IF ! time - t_bounce > t_bounce_lagr_chg
END IF ! lagr == 'ye'

!-----------------------------------------------------------------------
!  Check criteria for switching m_grid to 'no'
!-----------------------------------------------------------------------

IF ( m_grid == 'ye' ) THEN
  IF ( t_bounce > zero  .and.  time - t_bounce > t_bounce_mgrd_chg ) THEN
    m_grid                  = 'no'
    IF ( myid == 0 ) WRITE (nprint,103)
    WRITE (nlog,103)
  END IF ! time - t_bounce > t_bounce_mgrd_chg
END IF ! lagr == 'ye'

!-----------------------------------------------------------------------
!  Check criteria for switching time step criteria
!-----------------------------------------------------------------------

IF ( t_bounce > zero  .and.  time - t_bounce > 0.005 ) THEN
  tcntrl(21)                = 5.d-02
  tcntrl(22)                = 5.d-02
  tcntrl(23)                = 5.d-02
  tcntrl(24)                = 5.d-02
END IF ! t_bounce > zero  .and.  time - t_bounce > 0.005

!-----------------------------------------------------------------------
!  Check criteria for switching intnu_trns
!-----------------------------------------------------------------------

IF ( t_bounce > zero  .and.  time - t_bounce > 0.001 ) THEN
  intnu_trns                = 10
END IF

!-----------------------------------------------------------------------
!  Check criteria for switching on neutrino equilibration
!-----------------------------------------------------------------------

IF ( t_bounce > zero  .and.  time - t_bounce > t_equilibrate ) THEN
  nu_equil                  = 'ye'
END IF

!-----------------------------------------------------------------------
!  Check criteria for switching edit intervals
!-----------------------------------------------------------------------

IF ( ncycle == 10001 ) THEN

  intedc(1)                 = 200
  nedc(1)                   = 0
  nedc_r(1,:,:)             = 0

  intede(1)                 = 10000
  nede(1)                   = 0
  nede_r(1,:,:)             = 0

  DO i = 1,3
    intdma(i)               = 10000
    nedma(i)                = 0
    nedma_r(i,:,:)          = 0
  END DO

  DO i = 1,9
    intedh(i)               = 10000
    nedh(i)                 = 0
    nedh_r(i,:,:)           = 0
  END DO

  intdps(1)                 = 10000
  nedps(1)                  = 0
  nedps_r(1,:,:)            = 0

  DO i = 1,14
    intedu(i)               = 10000
    nedu(i)                 = 0
    nedu_r(i,:,:)           = 0
  END DO

  DO i = 1,6
    intedy(i)               = 10000
    nedy(i)                 = 0
    nedy_r(i,:,:)           = 0
  END DO

  DO i = 1,2
    intdsc(i)               = 10000
    nedsc(i)                = 0
    nedsc_r(i,:,:)          = 0
  END DO

  DO n = 1,nnu
    intedn(n)               = 10000
    nedn(n)                 = 0
    nedn_r(i,:,:)           = 0
  END DO

  DO i = 1,18
    intdng(i,:)             = 10000
    nedng(i,:)              = 0
    nedng_r(i,:,:,:)        = 0
  END DO ! 

  intprt                    = 1000
  nprt                      = 0

END IF ! ncycle = 20001

RETURN

END SUBROUTINE parameter_reset

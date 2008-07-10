SUBROUTINE e_advct_time_step( jr_min, jr_max )
!-----------------------------------------------------------------------
!
!    File:         e_advct_time_step
!    Module:       e_advct_time_step
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To select the time step restricted by neutrino energy advection.
!
!    Subprograms called:
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  e_advct_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY: nnu
USE numerical_module, ONLY: zero, epsilon

USE e_advct_module, ONLY : dtnph, t_cntrl_e_advct, dpsivmx, dt_nu_e_advct
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Inout variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min            ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max            ! maximum radial zone index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: n               ! neutrino flavor index

REAL(KIND=double), PARAMETER     :: dtmax = 1.d+20  ! times step without criteria
REAL(KIND=double), PARAMETER     :: dtmaxi = 1.d0/dtmax

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!  Time step controls
!-----------------------------------------------------------------------
!  Psi0 fraction change time step ( dpsi0nmax(n) ) due to advection of
!   neutrinos in energy.
!
!  t_cntrl_e_advct(n) is the maximum permitted abs( d(psi0)/( psi0 + psivmin )
!   due to advection of neutrinos in energy.
!  Psi0 fraction change time step control is bypassed if nnugp(n) = 0.
!  Psi0 fraction change time step control is bypassed if t_cntrl_e_advct(n) < 0.
!-----------------------------------------------------------------------

DO n = 1,nnu

!........Bypass if tcntrl(40+n) < 0, or no n-neutrinos, or ivc_x = 0

  IF ( nnugp(n) == 0  .or.  t_cntrl_e_advct(n) <= zero ) CYCLE

!-----------------------------------------------------------------------
!  Choose dt so that max relative change of psi0 not exceeded
!  dtmaxi sets dt_nu_e_advct(n) = 1.e20 if dmax = 0.0
!-----------------------------------------------------------------------

  dt_nu_e_advct(n) = t_cntrl_e_advct(n) * dtnph/( dpsivmx(n) + dtmaxi * t_cntrl_e_advct(n) * dtnph + epsilon )

END DO

RETURN
END SUBROUTINE e_advct_time_step

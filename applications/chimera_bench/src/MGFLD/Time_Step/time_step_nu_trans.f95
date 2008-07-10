SUBROUTINE time_step_nu_trans( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, dt, jdt, dtime_nutrns, nnu )
!-----------------------------------------------------------------------
!
!    File:         time_step_nu_trans
!    Module:       time_step_nu_trans
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To select the time step restricted by neutrino source and
!       transport when hydro subcycling is implemented.
!
!    Subprograms called:
!
!    Input arguments:
!  jr_min       : minimum radial index
!  jr_max       : maximum radial index
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  dtime_matter : minimum times step by hydro and other non-neutrino processes
!  nnu          : neutrino flavor array extent
!
!    Output arguments:
!  dt           : minimum transport time step restrictions for criteria
!                  i, radial ray ij_ray, ik_ray (s)
!  jdt          : radial zone causing minimum time step for criteria i,
!                  radial ray ij_ray, ik_ray
!  dtime_nutrns : minimum times step by neutrino source and transport
!                  criteria
!
!    Include files:
!  kind_module, numerical_module
!  cycle_module, incrmnt_module, mdl_cnfg.cmn, nu_dist.cmn,
!  nu_energy_grid_module, prb_cntl.cmn, t_cntrl.cmn
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: zero, epsilon

USE cycle_module, ONLY: nutrans_trns
USE incrmnt_module, ONLY: dtmpnn, dye
USE mdl_cnfg_module, ONLY: dr, rho, t, ye
USE nu_dist_module, ONLY: psi0
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE t_cntrl_module, ONLY: dtnmhn_trans, dtnph_trans, &
& dpsisp, psimin, rdtmax, tcntrl, delt_trans

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min           ! minimum radial index
INTEGER, INTENT(in)              :: jr_max           ! maximum radial index
INTEGER, INTENT(in)              :: ij_ray           ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray           ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nnu              ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)           :: jdt ! zone causing dt

REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim) :: dt  ! minimum allowed time step
REAL(KIND=double), INTENT(out)   :: dtime_nutrns    ! minimum times step by neutrino source and transport criteria

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i                ! do index
INTEGER                          :: j                ! radial zone index
INTEGER                          :: n                ! neutrino flavor index

REAL(KIND=double), PARAMETER     :: dtmax = 1.d+20   ! times step without criteria
REAL(KIND=double), PARAMETER     :: dtmaxi = 1.d-20  ! dtmax^-2
REAL(KIND=double), PARAMETER     :: dtn_mult = 1.d+1 ! elapsed time cannot exceed dtn_mult*dtnmh
REAL(KIND=double)                :: d                ! working variable
REAL(KIND=double)                :: dmax             ! maximum relative change in a quantity

REAL(KIND=double)                :: dtnph_t          ! time step restricted by source and transport processes (working variable)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nutrans_trns = false
!
!  If nutrans_trns = false the source and transport step were not
!   implemented during the current cycle, and new source and transport
!   time steps are therefore not computed.
!-----------------------------------------------------------------------

IF ( .not. nutrans_trns ) RETURN

!-----------------------------------------------------------------------
!
!                      \\\\\ INITIALIZE /////
!
!  Initialize ncynu_trans, dtnmhn_trans, dtnph_trans, and set nutrans_trns
!   to false since a source and transport step has been been taken this
!   cycle. The above parameters will be reset below and in subroutine
!   time_step_select.
!-----------------------------------------------------------------------

dtnmhn_trans          = dtnph_trans
dtime_nutrns          = dtmax

!----------------------------------------------------------------------!
!
!             \\\\\ HYDRO AND TRANSPORT CYCLING /////
!
!
!                  intnu_trns > 1  ... then
!
!  ith neutrino-transport is implemented at least once every intnu_trns
!   hydro updates, more often if intnu_trns hydro time steps exceed the
!   maximum allowed neutrino time step.
!
!  ncynu_trns = cycle counter for the hydro subcycling within a neutrino
!   transport step.
!
!  ncynu_trns = 1: ncynu_trns = intnu_trns during the preceding cycle
!   and has therefore been reset to 1. Neutrino transport has therefore
!   been implemented during the preceding cycle. Time steps for neutrino
!   transport are updated, and their values divided by intnu are stored
!   in array dts to be used as time step criteria for hydro until neutrino
!   transport is again implemented and new time step criteria determined.
!
!  ncynu_trns > 1: ncynu_trns < intnu_trns during the preceding cycle and
!   has therefore not been reset to 1. The neutrino transport has therefore
!   not been implemented during the preceding cycle. The neutrino transport
!   time step criteria stored in array dts are used  to limit the hydro time
!   steps.
!
!                  intnu_trns < 0  ... then
!
!  Separate timesteps dtime_matter and dtnph_trans are computed for hydro
!   and neutrino transport, respectively.
!
!  dtnph_trans < dtime_matter : dtime_matter is set equal to dtnph_trans
!   and neutrino transport is implemented.
!
!  dtnph_trans > dtime_matter : neutrino transport is bypassed until the
!   accumulated time since the last neutrino transport step could exceed
!   dtnph_trans  by a specified amount in the subsequent time step. Thus,
!   if delt_i is the accumulated time from the last implementation of
!   neutrino transport to the beginning of the current time cycle, neutrino
!   transport is implemented in the current time cycle with dtnph_trans
!   set equal to the accumulated time if
!
!               delt_i + rdtmax*dtnph_t > dtnph_trans .
!
!  To ensure that the accumulated time will not exceed the predicted ith
!   neutrino process time step dtnph_trans, rdtmax must be given by
!
!               rdtmax = 1 + tcntrl(10) .
!
!  Setting rdtmax to a negative number in the data file will instruct the
!   code to reset rdtmax to the above value.
!
!.......................................................................
!
!        Timestep criteria.
!
!  dt(i,ij_ray,ik_ray)    : minimum time step given by criterion i for
!                    radial ray ij_ray, ik_ray
!  jdt(i,ij_ray,ik_ray)   : zone at which minimum timestep given by
!                    criterion i occurs for radial ray ij_ray, ik_ray
!
!  dt(11,ij_ray,ik_ray)   : Net temperature change time step due to
!                    source and transport; tcntrl(11) is the maximum
!                    permitted abs( dtmpnn(j,n,ij_ray,ik_ray)/t(j) ).
!
!  dt(16,ij_ray,ik_ray)   : Net electron fraction change time step due
!                    source and transport; tcntrl(16) is the maximum
!                    permitted abs( dye(j,n,ij_ray,ik_ray)/ye(j) ).
!
!  dt(20+n,ij_ray,ik_ray) : Psi0 fraction change time step due to source
!                    and transport of n-neutrinos; tcntrl(20+n) is the
!                    maximum permitted
!                    abs( d(psi0(:,:,n)/( psi0(:.:.n) + psimin ) ).
!        
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize time step restrictions; jdt(i,ij_ray,ik_ray), i = 21 - 24
!   are set in subroutine nu_trans.
!-----------------------------------------------------------------------

DO i = 11,24
  dt(i,ij_ray,ik_ray)       = dtmax
END DO

DO i = 11,20
  jdt(i,ij_ray,ik_ray)      = 0
END DO

!-----------------------------------------------------------------------
!  If jdt(50,ij_ray,ik_ray) = -1, time step is set equal to tcntrl(50);
!   otherwise tcntrl(50) is the maximum allowable time step.
!-----------------------------------------------------------------------

IF ( jdt(50,ij_ray,ik_ray) == -1 ) THEN
  dtnph_t                   = tcntrl(50)
  dtime_nutrns              = tcntrl(50)
  RETURN
END IF

!-----------------------------------------------------------------------
!
!            \\\\\ NET TEMPERATURE CHANGE CRITERIA /////
!
!  Net temperature change time step control due to source and transport
!
!               ***( dt(11,ij_ray,ik_ray) )***
!
!  Net temperature change due to source and transport is stored in
!   dtmpnn(j,1,ij_ray,ik_ray).
!  tcntrl(11) is the maximum permitted abs( dtmpnn(j,1,ij_ray,ik_ray)/t(j) ).
!  Net temperature change time step control is bypassed if tcntrl(10+n) < 0.
!  Net temperature change time step control is bypassed if nnugpmx = 0.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine restriction if tcntrl(11) > 0 and neutrinos presemt
!-----------------------------------------------------------------------

IF ( tcntrl(11) > zero  .and.  nnugpmx > 0 ) THEN

!-----------------------------------------------------------------------
!  Find maximum relative change of temperature
!-----------------------------------------------------------------------

  dmax                      = zero

  DO j = jr_min,jr_max
    d                       = DABS( dtmpnn(j,1,ij_ray,ik_ray)/t(j) )
    IF ( d >= dmax ) THEN
      dmax                  = d
      jdt(11,ij_ray,ik_ray) = j
    END IF ! d > dmax
  END DO

!-----------------------------------------------------------------------
!  Choose dt so that max change restriction in t is not exceeded 
!     dtmaxi sets dt(11,ij_ray,ik_ray) = 1.e20 if dmax = 0.0
!-----------------------------------------------------------------------

  dt(11,ij_ray,ik_ray)      = tcntrl(11) * dtnmhn_trans                &
&                           / ( dmax + dtmaxi * tcntrl(11) * dtnmhn_trans + epsilon )

END IF ! tcntrl(11) > zero  .and.  nnugpmx > 0

!-----------------------------------------------------------------------
!
!          \\\\\ NET ELECTRON FRACTION CHANGE CRITERIA /////
!
!  Net electron fraction change time step due to source and transport,
!   of n-neutrinos.
!
!                   ***( dt(16,ij_ray,ik_ray) )***
!
!  tcntrl(16) is the maximum permitted abs( dye(j,n,ij_ray,ik_ray)/ye(j) )
!  Net electron fraction change time step control is bypassed if
!   tcntrl(16) < 0.
!  Net electron fraction change time step control is bypassed if
!   nnugpmx = 0.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine restriction if tcntrl(16) > 0 and neutrinos presemt
!-----------------------------------------------------------------------

IF ( tcntrl(16) > zero  .and.  nnugpmx > 0 ) THEN

!-----------------------------------------------------------------------
!  Get maximum relative change of electron fraction
!-----------------------------------------------------------------------

  dmax                      = zero
  DO j = jr_min,jr_max
    d                       = DABS( dye(j,1,ij_ray,ik_ray)/ye(j) )
    IF ( d > dmax ) THEN
      dmax                  = d
      jdt(16,ij_ray,ik_ray) = j
    END IF ! d > dmax
  END DO

!-----------------------------------------------------------------------
!  Choose dt so that max change in ye not exceeded
!     dtmaxi sets dt(16,ij_ray,ik_ray) = 1.e20 if dmax = 0.0
!-----------------------------------------------------------------------

  dt(16,ij_ray,ik_ray)      = tcntrl(16) * dtnmhn_trans                &
&                           / ( dmax + dtmaxi * tcntrl(16) * dtnmhn_trans + epsilon )

END IF ! tcntrl(16) > zero  .and.  nnugpmx > 0

!-----------------------------------------------------------------------
!
!               \\\\\ NET PSI0 CHANGE CRITERIA /////
!
!  Net psi0 fraction change time step due to source and transport.
!
!               ***( dt(20+n,ij_ray,ik_ray) )***
!
!  tcntrl(20+n) is the maximum permitted abs( d(psi0)/psi0 ) due to
!   emsision, absorption, scattering, production and transport.
!  Psi0 fraction change time step control is bypassed if tcntrl(20+n) < 0.
!  Psi0 fraction change time step control is bypassed if nnugpmx = 0.
!-----------------------------------------------------------------------

DO n = 1,nnu

!-----------------------------------------------------------------------
!  Bypass if tcntrl(20+n) < 0 or no n-neutrinos
!-----------------------------------------------------------------------

  IF ( nnugpmx == 0  .or.  tcntrl(20+n) <= zero ) CYCLE

!-----------------------------------------------------------------------
!  Time step due to maximum relative change of psi0
!  Choose dt so that max change in psi0 not exceeded
!  dtmaxi sets dt(20+n,ij_ray,ik_ray) = 1.e20 if dmax = 0.0
!-----------------------------------------------------------------------

  dt(20+n,ij_ray,ik_ray) = tcntrl(20+n) * dtnmhn_trans                 &
&                        / ( dpsisp(n) + dtmaxi * tcntrl(20+n) * dtnmhn_trans + epsilon )

END DO


!-----------------------------------------------------------------------
!
!               \\\\\ SELECT MINIMUM TIME STEP /////
!
!-----------------------------------------------------------------------

DO i = 11,24
  dtime_nutrns           = DMIN1( dtime_nutrns, dt(i,ij_ray,ik_ray) )
END DO

RETURN
END SUBROUTINE time_step_nu_trans

!-----------------------------------------------------------------------
!    Module:       cycle_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE cycle_module

USE kind_module
SAVE

!-----------------------------------------------------------------------
!  Computational cycles.
!-----------------------------------------------------------------------
!  ncycle : the cycle number of the calculation.
!
!  ncymax : probelm termination crierion. Calculation terminated when
!   ncycle > ncymax.
!
!  nrst   : the cycle number at problem initiation.
!-----------------------------------------------------------------------

INTEGER                                        :: ncycle
INTEGER                                        :: ncymax
INTEGER                                        :: nrst

!-----------------------------------------------------------------------
!  Hydro subcycling
!-----------------------------------------------------------------------
! iintnu_trns  :  hydro subcycling parameters
!
!   > 0: intnu_i is the number of 'hydro' cycles per neutrino cycle at
!    the current time.
!   < 0: separate timesteps dtnph and dtnphn_i are computed for 'hydro'
!    and neutrino step, respectively.
!
!   dtnphn_i   < dtnph: dtnph is set equal to dtnphn_i and neutrino step is
!
!   dtnpnn_i > dtnph: neutrino process i is bypassed until the accumulated
!    time, delt_i, since implemented. If the last step plus the current 
!    'hydro' time step could exceed dtnphn_i in the next cycle, i.e.,
!    untildtnphn_i. When the above delt_i + dtnph*rdtmax > inequality is
!    satisfied, neutrino transport step is called with dtnphn_i = delt_i.
!
!  intnur_trns : hydro subcycling parameters for the preceding time step.
!   i.e., the number of hydro updates per neutrino step at the last
!   neutrino transport step.
!
!  ncynu_trns  : the number of hydro updates since the last neutrino step
!-----------------------------------------------------------------------

LOGICAL                                        :: nutrans_trns

INTEGER                                        :: intnu_trns
INTEGER                                        :: intnur_trns
INTEGER                                        :: ncynu_trns

END module cycle_module

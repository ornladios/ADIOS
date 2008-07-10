SUBROUTINE unpack_e_advct_keys( nnu, i_e_advct_data, d_e_advct_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_e_advct_keys
!    Module:       unpack_e_advct_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/26/04
!
!    Purpose:
!      To unpack the neutrino e_advection key arrays and restore the values
!       to the appropriate variables in e_advct_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nnu            : neutrino flavor dimension
!  i_e_advct_data : integer array of edit keys
!  d_e_advct_data : real*8 array of edit keys
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  e_advct_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE e_advct_module, ONLY : ivc_x, ivc_y, ivc_z, tau_advct, t_cntrl_e_advct, &
& psivmin, rhomin_y_eadvect, rhomin_z_eadvect

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nnu          ! neutrino flavor dimension

INTEGER, INTENT(in), DIMENSION(5)                 :: i_e_advct_data  ! integer array of e_advect keys

REAL(KIND=double), INTENT(in), DIMENSION(5+2*nnu) :: d_e_advct_data  ! 64 bit real array of e_advect keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: n             ! neutrino flavor index

!-----------------------------------------------------------------------
!
!                \\\\\ UNPACK E_ADVECTION KEYS /////
!
!-----------------------------------------------------------------------


ivc_x                    = i_e_advct_data(1)
ivc_y                    = i_e_advct_data(2)
ivc_z                    = i_e_advct_data(3)

tau_advct                = d_e_advct_data(1)
rhomin_y_eadvect         = d_e_advct_data(2)
rhomin_z_eadvect         = d_e_advct_data(3)

DO n = 1,nnu
  t_cntrl_e_advct(n)     = d_e_advct_data(5+0*nnu+n)
  psivmin(n)             = d_e_advct_data(5+1*nnu+n)
END DO

RETURN
END SUBROUTINE unpack_e_advct_keys

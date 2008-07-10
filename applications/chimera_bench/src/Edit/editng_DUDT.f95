SUBROUTINE editng_DUDT( n, jr_min, jr_max )
!-----------------------------------------------------------------------
!
!    File:         editng_DUDT
!    Module:       editng_DUDT
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/17/05
!
!    Purpose:
!      To edit integral n-type neutrino transport data pertaining to
!       neutrino-nucleon inelastic scattering (NNS).
!
!    Subprograms call:
!  sctedns    : Computes quantities needed for editing NNS rates
!
!    Input arguments:
!  n          : neutrino flavor
!  jr_min     : inner radial zone of region for which configuration edit is to be made.
!  jr_max     : outer radial zone of region for which configuration edit is to be made.
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!      none
!
!    Include files:
!      kind_module, array_module, numerical_module, physcnst_module
!      edit_module, mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : half

USE edit_module, ONLY : nprint, head, dudt_ABEM, dudt_Brem, dudt_NAS, &
& dudt_NES, dudt_NNS, dudt_PR, dudt_NET
USE mdl_cnfg_module, ONLY : rho, rstmss
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: n             ! neutrino flavor
INTEGER, INTENT(in)               :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)               :: jr_max        ! maximum radial zone index

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: j             ! radial zone index

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1x)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*')/)
 1501 FORMAT (30x,'Neutrino-Nucleon Net Energy-Matter Transfer Rates (MeV/nucleon)'/)
 1503 FORMAT ('   j     rho       rstmss   e_nu NET   e_nu ABEM eb_nu ABEM &
 & e_nu NES   eb_nu NES  x_nu NES   xb_nu NES  e_nu Pair  eb_nu Pair x_nu Pair &
 & xb_nu Pair e_nu NNS   eb_nu NNS  x_nu NNS   xb_nu NNS  e_nu BRM   eb_nu BRM &
 & x_nu BRM   xb_nu BRM  e_nu NAS   eb_nu NAS  x_nu NAS   xb_nu NAS'/)
 1505 FORMAT (1x,i4,25(es11.3))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) RETURN

!-----------------------------------------------------------------------
!  Compute dudt_NET, the net energy transferred to matter by neutrinos
!   (Mev/nucleon)
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  dudt_NET(j)    = dudt_ABEM(j,1) + dudt_ABEM(j,2)                                                 &
&                + dudt_NES(j,1)  + dudt_NES(j,2)  + 2.d0 * dudt_NES(j,3)  + 2.d0 * dudt_NES(j,4)  &
&                + dudt_PR(j,1)   + dudt_PR(j,2)   + 2.d0 * dudt_PR(j,3)   + 2.d0 * dudt_PR(j,4)   &
&                + dudt_NNS(j,1)  + dudt_NNS(j,2)  + 2.d0 * dudt_NNS(j,3)  + 2.d0 * dudt_NNS(j,4)  &
&                + dudt_Brem(j,1) + dudt_Brem(j,2) + 2.d0 * dudt_Brem(j,3) + 2.d0 * dudt_Brem(j,4) &
&                + dudt_NAS(j,1)  + dudt_NAS(j,2)  + 2.d0 * dudt_NAS(j,3)  + 2.d0 * dudt_NAS(j,4)
END DO

!-----------------------------------------------------------------------
!
!                         \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

WRITE (nprint,1)
WRITE (nprint,3) head
WRITE (nprint,5)
CALL date_and_time_print( nprint )
WRITE (nprint,1501)
WRITE (nprint,1503)

DO j = jr_max,jr_min,-1
  WRITE (nprint,1505) j,rho(j), half*(rstmss(j)+rstmss(j-1)), dudt_NET(j), &
& dudt_ABEM(j,1), dudt_ABEM(j,2), dudt_NES(j,1), dudt_NES(j,2), dudt_NES(j,3), &
& dudt_NES(j,4), dudt_PR(j,1), dudt_PR(j,2), dudt_PR(j,3), dudt_PR(j,4), &
& dudt_NNS(j,1), dudt_NNS(j,2), dudt_NNS(j,3), dudt_NNS(j,4), dudt_Brem(j,1), &
& dudt_Brem(j,2), dudt_Brem(j,3), dudt_Brem(j,4), dudt_NAS(j,1), dudt_NAS(j,2), &
& dudt_NAS(j,3), dudt_NAS(j,4)
END DO ! j

RETURN
END SUBROUTINE editng_DUDT

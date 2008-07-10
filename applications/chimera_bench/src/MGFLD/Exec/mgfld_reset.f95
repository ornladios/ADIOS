SUBROUTINE mgfld_reset( ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         mgfld_reset
!    Module:       mgfld_reset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/09/09
!
!    Purpose:
!      To reset the neutrino rate tables.
!
!    Subprograms called:
!  abemset   : recomputes, if necessary, table of neutrino emission and
!               absorption rates
!  bremset   : recomputes, if necessary, table of nucleon-nucleus
!               bremsstrahlung rates
!  pairset   : recomputes, if necessary, table of pair annihilation rates
!  scataset  : recomputes, if necessary, table of neutrino-nucleus
!               inelastic scattering rates
!  scateset  : recomputes, if necessary, table of neutrino-electron
!               elastic scattering rates
!  scatiset  : recomputes, if necessary, table of neutrino-nucleon and
!               nucleus isoenergetic scattering rates
!  scatnnset : recomputes, if necessary, table of neutrino-nucleon
!               inelastic scattering rates
!  scatnAset : recomputes, if necessary, table of neutrino-nucleus
!               inelastic scattering rates
!  scatnset  : recomputes, if necessary, table of neutrino-nucleon
!               elastic scattering rates
!
!    Input arguments:
!  ij_ray    : j-index of a radial ray
!  ik_ray    : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  mdl_cnfg_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double

USE mdl_cnfg_module, ONLY: jr_min, jr_max, rho, t, ye
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: j             ! zone index
INTEGER                          :: jr_maxp       ! jr_max + 1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ RESET NEUTRINO RATE TABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Recompute neutrino rates if nnugpmx > 0
!-----------------------------------------------------------------------

IF ( nnugpmx > 0 ) THEN

  jr_maxp          = jr_max + 1

  DO j = jr_min,jr_maxp
    CALL abemset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
    CALL bremset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
    CALL pairset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
    CALL scataset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
    CALL scateset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
    CALL scatiset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
    CALL scatnAset( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
    CALL scatnset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
    CALL scatnnset( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO

END IF !  nnugpmx > 0

RETURN
END SUBROUTINE mgfld_reset

SUBROUTINE enua_cal( jr_min, jr_max )
!-----------------------------------------------------------------------
!
!    File:         enua_cal
!    Module:       enua_cal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/00
!
!    Purpose:
!      To compute e_nu/a.
!
!    Variables that must be passed through common:
!  agr_nu(j)     : GR lapse function
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min        : inner zone for which calculation w is to be made
!  jr_max        : outer zone for which calculation w is to be made
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  nnugp(n)      : number of neutrino energy groups for n-neutrinos
!  unubi(k)      : energy of outer edge of neutrino energy zone k as
!                   measured by an observer at infinity
!  unui(k)       : zone centered energy of neutrino energy zone k as
!                   measured by an observer at infinity
!  dunui(i)      : unubi(k+1) - unubi(k)
!  agra_nu(j)    : lapse function at zone j at time m
!  agrajmh_nu(j) : lapse function at zone j-1/2 at time m
!
!    Output arguments (common):
!  unuea(j,k)    : zone centered energy of neutrino energy zone k at
!                   radial zone j at time m
!  unubea(j,k)   : energy of outer edge of neutrino energy zone k at
!                   radial zone j at time m
!  dunuea(j,k)   : unube(j,k+1) - unube(j,k)
!  unua(j,k)     : zone centered energy of neutrino energy zone k at
!                   radial zone j-1/2 at time m
!  unuba(j,k)    : energy of outer edge of neutrino energy zone k at
!                   radial zone j-1/2 at time m
!  dunua(j,k)    : unub(j,k+1) - unub(j,k)
!  ncoefaa(j,k)  : number of neutrino states per unit volume at radial
!                   zone j energy unu(j,k)
!  ecoefaa(j,k)  : energy sum neutrino states per unit volume at radial
!                   zone j-1/2 energy unu(j,k)
!  ecoefaea(j,k) : energy sum neutrino states per unit volume at radial
!                   zone j energy unub(j,k)
!
!    Include files:
!  numerical_module
!  nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : ncoef, ecoef

USE nu_dist_module, ONLY : unuea, dunuea, unubea, unua, dunua, unuba, &
& agra_nu, agrajmh_nu, ncoefaa, ecoefaa, ecoefaea
USE nu_energy_grid_module, ONLY : nnugp, unui, unubi, dunui, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(IN)              :: jr_max        ! maximum radial zone index

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: kp            ! nnugpmx+1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Incorporate GR corrections to neutrino energies
!-----------------------------------------------------------------------

kp                               = nnugpmx + 1

DO k = 1,nnugpmx

  unuea(jr_min:jr_max,k)         = unui(k) /agra_nu(jr_min:jr_max)
  dunuea(jr_min:jr_max,k)        = dunui(k)/agra_nu(jr_min:jr_max)
  unubea(jr_min:jr_max,k)        = unubi(k)/agra_nu(jr_min:jr_max)

  unua(jr_min:jr_max,k)          = unui(k) /agrajmh_nu(jr_min:jr_max)
  dunua(jr_min:jr_max,k)         = dunui(k)/agrajmh_nu(jr_min:jr_max)
  unuba(jr_min:jr_max,k)         = unubi(k)/agrajmh_nu(jr_min:jr_max)

END DO ! k

unubea(jr_min:jr_max,kp)          = unubi(kp)/agra_nu(jr_min:jr_max)
unuba(jr_min:jr_max,kp)           = unubi(kp)/agrajmh_nu(jr_min:jr_max)

ncoefaa (jr_min:jr_max,1:nnugpmx) = ncoef * unua (jr_min:jr_max,1:nnugpmx)**2 * dunua (jr_min:jr_max,1:nnugpmx)
ecoefaa (jr_min:jr_max,1:nnugpmx) = ecoef * unua (jr_min:jr_max,1:nnugpmx)**3 * dunua (jr_min:jr_max,1:nnugpmx)
ecoefaea(jr_min:jr_max,1:nnugpmx) = ecoef * unuea(jr_min:jr_max,1:nnugpmx)**3 * dunuea(jr_min:jr_max,1:nnugpmx)

ncoefaa (jr_max+1,1:nnugpmx)      = ncoefaa (jr_max,1:nnugpmx)
ecoefaa (jr_max+1,1:nnugpmx)      = ecoefaa (jr_max,1:nnugpmx)
ecoefaea(jr_max+1,1:nnugpmx)      = ecoefaea(jr_max,1:nnugpmx)
ecoefaea(1       ,1:nnugpmx)      = ecoefaea(2     ,1:nnugpmx)


RETURN
END SUBROUTINE enua_cal

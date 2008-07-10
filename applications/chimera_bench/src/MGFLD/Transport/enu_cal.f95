SUBROUTINE enu_cal( jr_min, jr_max )
!-----------------------------------------------------------------------
!
!    File:         enu_cal
!    Module:       enu_cal
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
!  agr_nu(j)   : GR lapse function
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min       : inner zone for which calculation w is to be made
!  jr_max       : outer zone for which calculation w is to be made
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  nnugp(n)     : number of neutrino energy groups for n-neutrinos
!  unubi(k)     : energy of outer edge of neutrino energy zone k as
!                  measured by an observer at infinity
!  unui(k)      : zone centered energy of neutrino energy zone k as
!                  measured by an observer at infinity
!  dunui(i)     : unubi(k+1) - unubi(k)
!  agr_nu(j)    : lapse function at zone j at time m
!  agrjmh_nu(j) : lapse function at zone j-1/2 at time m
!
!    Output arguments (common):
!  unue(j,k)    : zone centered energy of neutrino energy zone k at
!                  radial zone j at time m
!  unube(j,k)   : energy of outer edge of neutrino energy zone k at
!                  radial zone j at time m
!  dunue(j,k)   : unube(j,k+1) - unube(j,k)
!  unu(j,k)     : zone centered energy of neutrino energy zone k at
!                  radial zone j-1/2 at time m
!  unub(j,k)    : energy of outer edge of neutrino energy zone k at
!                  radial zone j-1/2 at time m
!  dunu(j,k)    : unub(j,k+1) - unub(j,k)
!  ncoefa(j,k)  : number of neutrino states per unit volume at radial
!                  zone j energy unu(j,k)
!  ecoefa(j,k)  : energy sum neutrino states per unit volume at radial
!                  zone j-1/2 energy unu(j,k)
!  ecoefae(j,k) : energy sum neutrino states per unit volume at radial
!                  zone j energy unub(j,k)
!
!    Include files:
!  numerical_module
!  nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : ncoef, ecoef

USE nu_dist_module, ONLY : unue, dunue, unube, unu, dunu, unub, &
& agr_nu, agrjmh_nu, ncoefa, ecoefa, ecoefae
USE nu_energy_grid_module, ONLY : nnugp, unui, unubi, dunui, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index

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
  
  unue (jr_min:jr_max,k)         = unui(k) /agr_nu(jr_min:jr_max)
  dunue(jr_min:jr_max,k)         = dunui(k)/agr_nu(jr_min:jr_max)
  unube(jr_min:jr_max,k)         = unubi(k)/agr_nu(jr_min:jr_max)

  unu (jr_min:jr_max,k)          = unui(k) /agrjmh_nu(jr_min:jr_max)
  dunu(jr_min:jr_max,k)          = dunui(k)/agrjmh_nu(jr_min:jr_max)
  unub(jr_min:jr_max,k)          = unubi(k)/agrjmh_nu(jr_min:jr_max)

END DO ! k

unube(jr_min:jr_max,kp)          = unubi(kp)/agr_nu   (jr_min:jr_max)
unub (jr_min:jr_max,kp)          = unubi(kp)/agrjmh_nu(jr_min:jr_max)

ncoefa (jr_min:jr_max,1:nnugpmx) = ncoef * unu (jr_min:jr_max,1:nnugpmx)**2 * dunu (jr_min:jr_max,1:nnugpmx)
ecoefa (jr_min:jr_max,1:nnugpmx) = ecoef * unu (jr_min:jr_max,1:nnugpmx)**3 * dunu (jr_min:jr_max,1:nnugpmx)
ecoefae(jr_min:jr_max,1:nnugpmx) = ecoef * unue(jr_min:jr_max,1:nnugpmx)**3 * dunue(jr_min:jr_max,1:nnugpmx)

ncoefa (jr_max+1,1:nnugpmx)      = ncoefa (jr_max,1:nnugpmx)
ecoefa (jr_max+1,1:nnugpmx)      = ecoefa (jr_max,1:nnugpmx)
ecoefae(jr_max+1,1:nnugpmx)      = ecoefae(jr_max,1:nnugpmx)
ecoefae(1       ,1:nnugpmx)      = ecoefae(2     ,1:nnugpmx)

RETURN
END SUBROUTINE enu_cal

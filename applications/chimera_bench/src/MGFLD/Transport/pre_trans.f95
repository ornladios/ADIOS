SUBROUTINE pre_trans( jr_min, jr_max, rho, r, nx, nnu )
!-----------------------------------------------------------------------
!
!    File:         pre_trans
!    Module:       pre_trans
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/10/00
!
!    Purpose:
!      To compute quantities needed in the neutrino transport
!       modules.
!
!    Subprograms called:
! gamgr_nu_cal : update the GR gamma's
!
!    Input arguments:
!
!  jr_min       : inner radial zone number
!  jr_max       : outer radial zone number
!  rho          : density (g cm^{-3}
!  r            : radius (cm)
!  nx           : -array extent
!  nnu          : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!
!  nnugp(n)     : number of energy zones for neutrinos of type n
!  unu(j,k)     : midpoint energy of energy zone k at radial zone j
!  dunu(j,k)    : energy width of energy zone k at radial zone j
!
!    Output arguments (common):
!
!  ncoefa(j,k)  : 4.*pi/((h*!)**3)*w**2*dw
!  ecoefa(j,k)  : 4.*pi*ergmev/((h*!)**3)*w**3*dw at j-1/2
!  ecoefae(j,k) : 4.*pi*ergmev/((h*!)**3)*w**3*dw at j
!  area(j)      : area of radial zone j
!  vol(j)       : volume enclosed by radial zone j and j-1
!  drjmh(j)     : proper distance between r(j) and r(j-1)
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, half, one, frpi, ncoef, ecoef
USE physcnst_module, ONLY : rmu

USE mdl_cnfg_module, ONLY : dmrst
USE nu_dist_module, ONLY : stwt, unu, dunu, unue, dunue, ncoefa, ecoefa, ecoefae,  &
& rjmh, gamgr_nu, area, areajmh, vol, drjmh, drjmh_inv
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: nx            ! radial array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)         :: rho      ! density (cm^{-3})

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx)      :: r        ! radial zone radii (cm)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Functions of the neutrino energy
!-----------------------------------------------------------------------

ncoefa (jr_min:jr_max,:) = ncoef * unu (jr_min:jr_max,:)**2 * dunu (jr_min:jr_max,:)
ecoefa (jr_min:jr_max,:) = ecoef * unu (jr_min:jr_max,:)**3 * dunu (jr_min:jr_max,:)
ecoefae(jr_min:jr_max,:) = ecoef * unue(jr_min:jr_max,:)**3 * dunue(jr_min:jr_max,:)

!-----------------------------------------------------------------------
!  Geometrical factors
!-----------------------------------------------------------------------

CALL gamgr_nu_cal( jr_min, jr_max) 

area(jr_min-1)           = zero
drjmh(jr_min-1)          = half * ( r(jr_min) )/gamgr_nu(jr_min)
r(jr_max+1)              = r(jr_max) + ( r(jr_max) - r(jr_max-1) )

rjmh     (jr_min:jr_max) = half * ( r(jr_min:jr_max) + r(jr_min-1:jr_max-1) )
area     (jr_min:jr_max) = frpi * r(jr_min:jr_max) * r(jr_min:jr_max)
areajmh  (jr_min:jr_max) = frpi * rjmh(jr_min:jr_max) * rjmh(jr_min:jr_max)
vol      (jr_min:jr_max) = dmrst(jr_min:jr_max)/rho(jr_min:jr_max)
drjmh    (jr_min:jr_max) = half * ( r(jr_min+1:jr_max+1) - r(jr_min-1:jr_max-1) )/gamgr_nu(jr_min:jr_max)
drjmh_inv(jr_min:jr_max) = one/drjmh(jr_min:jr_max)

rjmh(jr_max+1)           = r(jr_max) + half * ( r(jr_max) - r(jr_max-1) )

RETURN
END SUBROUTINE pre_trans

!-----------------------------------------------------------------------
!    Module:       it_tol_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE it_tol_module

USE kind_module
SAVE

!-----------------------------------------------------------------------
!  Neutrino transport tolerances
!-----------------------------------------------------------------------
!  iternu  : maximum number of iterations to attempt in order to obtain
!   convergence of the neutrino transport variables (i.e., t, ye, and the
!   psi0's). If iternu = 1, the variables are assumed to have converged
!   after the first iteration attempt.
!
!  itfail  : iteration failure switch.
!
!     itfail : 0 - iteration failure does not stop calculation.
!     itfail : 1 - iteration failure stops calculation.
!
!  tolnut  : Temperature convergence parameter for neutrino transport. The
!   criterion for temperature convergence is that
!
!     abs(dt/t) < tolnut .
!
!  tolnuye : Electron fraction convergence parameter for neutrino transport.
!   The criterion for electron fraction convergence is that
!
!     abs( dye/ye ) < tolnuye .
!
!  tolnupsi, tolpsimin : zero-monent neutrino occupation distribution
!   convergence parameters for neutrino transport. The criteria for neutrino
!   occupation convergence is that
!
!     abs( dpsi0/max( psi0, tolpsimin ) < tolnupsi .
!
!  a_prec  :  recompute neutrino rates if
!
!     dabs( agr(j) - agr_prev(j) )/agr(j) > a_prec
!-----------------------------------------------------------------------

INTEGER                                         :: iternu
INTEGER                                         :: itfail

REAL(KIND=double)                               :: tolnut
REAL(KIND=double)                               :: tolnuye
REAL(KIND=double)                               :: tolnupsi
REAL(KIND=double)                               :: tolpsimin
REAL(KIND=double)                               :: a_prec

!-----------------------------------------------------------------------
!  Hydrodynamic tolerances
!-----------------------------------------------------------------------
!  it_hy : maximum number of iterations to attempt in order to obtain a
!   convergent solution of the density, temperature, inverse radial metric
!   component, and gravitational mass due to hydrodynamics. If itnuc = 1,
!   the variables are assumed to have converged after the first iteration
!   attempt.
!
!  tolhy : convergence tolerance for hydrodynamic variables. The criteria
!   for convergence is that
!
!     abs(drho/rho)       < tolhy
!     abs(dt/t)           < tolhy
!     abs(dgamgr/gamgr)   < tolhy
!     abs(ddmgrv/dmgrv)   < tolhy
!
!     dabs( agr(j) - agr_prev(j) )/agr(j) > a_prec
!-----------------------------------------------------------------------

INTEGER                                         :: it_hy
REAL(KIND=double)                               :: tolhy

END module it_tol_module

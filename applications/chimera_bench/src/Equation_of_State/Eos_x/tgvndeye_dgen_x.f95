SUBROUTINE tgvndeye_dgen_x( j, ij_ray, ik_ray, rho, dr, e, ye, t_prev, t )
!-----------------------------------------------------------------------
!
!    File:         tgvndeye_dgen_x
!    Module:       tgvndeye_dgen_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/29/03
!
!    Purpose:
!        To compute the temperature given rho, s, and ye. Iteration
!         is by means of the bisection method if an initial guess
!         of the temperature is unavailable or is not within fraction
!         of the previous value, otherwise by Newton-Rhapson.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!  eqstt_x : interpolates eos quantities
!
!    Input arguments:
!  j       : radial zone index.
!  ij_ray  : index denoting the j-index of a specific radial ray
!  ik_ray  : index denoting the k-index of a specific radial ray
!  rho     : density (g/cm**3).
!  dr      : delta density (g/cm**3).
!  e       : inernal energy (ergs/g)
!  ye      : electron fraction.
!
!    Output arguments:
!  t       : temperature (K)
!  dt      : change in temperature (K)
!
!    Include files:
!  kind_module
!  eos_snc_x_module, edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE eos_snc_x_module, ONLY: aesvt
USE edit_module, ONLY: nread, nprint

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! index denoting the k-index of a specific radial ray

REAL(KIND=double), INTENT(in)    :: dr            ! delta density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: rho           ! density (g/cm**3)
REAL(KIND=double), INTENT(in)    :: e             ! internal energy (ergs/g)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction
REAL(KIND=double), INTENT(in)    :: t_prev        ! guess of the temperature

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: t             ! temperature (K)

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

REAL(KIND=double)                :: dum1          ! dummy variable
REAL(KIND=double)                :: dum2          ! dummy variable
REAL(KIND=double)                :: dum3          ! dummy variable
REAL(KIND=double)                :: dpdt          ! d(pressure)/dt
REAL(KIND=double)                :: dedt          ! d(energy)/dt

REAL(KIND=double)                :: dt            ! increment in temperature

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Temperature increment due to hydrodynamics
!-----------------------------------------------------------------------

 dt               = ( dr/( rho * ( rho + dr ) ) ) * t_prev              &
 &                * ( aesvt(j,1,ij_ray,ik_ray)/aesvt(j,2,ij_ray,ik_ray) )
 t                = dt + t_prev

 CALL eqstt_x( 1, j, ij_ray, ik_ray, rho + 0.5d0 * dr, 0.5d0 * ( t + t_prev ), &
 & ye, dum1, dum2, dpdt, dum3 )
 CALL eqstt_x( 2, j, ij_ray, ik_ray, rho + 0.5d0 * dr, 0.5d0 * ( t + t_prev ), &
 & ye, dum1, dum2, dedt, dum3)
 dt               = ( dr/( rho * ( rho + dr ) ) ) * 0.5d0 * ( t + t_prev ) * ( dpdt/dedt )
 t                = dt + t_prev

 CALL eqstt_x( 1, j, ij_ray, ik_ray,rho + 0.5d0 * dr, 0.5d0 * ( t + t_prev ), &
 & ye, dum1, dum2, dpdt, dum3 )
 CALL eqstt_x( 2, j, ij_ray, ik_ray,rho + 0.5d0 * dr, 0.5d0 * ( t + t_prev ), &
 & ye, dum1, dum2, dedt, dum3 )
 dt               = ( dr/( rho * ( rho + dr ) ) ) * 0.5d0 * ( t + t_prev ) * ( dpdt/dedt )
 t                = dt + t_prev

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE tgvndeye_dgen_x

SUBROUTINE tgvndeye_dgen_y( j, ji_ray, jk_ray, rho, dr, e, ye, t_prev, t )
!-----------------------------------------------------------------------
!
!    File:         tgvndeye_dgen_y
!    Module:       tgvndeye_dgen_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/19/05
!
!    Purpose:
!      To compute the temperature given rho, s, and ye. Iteration
!       is by means of the bisection method if an initial guess
!       of the temperature is unavailable or is not within fraction
!       of the previous value, otherwise by Newton-Rhapson.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!  eqstt_y : computes selected EOS quantities
!
!    Input arguments:
!  j       : y (angular) zone index.
!  ji_ray  : x (radial) index of a specific y (angular) ray
!  jk_ray  : z (azimuthal) index of a specific y (angular) ray
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
!  eos_snc_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE eos_snc_y_module, ONLY: aesvt

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! y (angular) zone index
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific y (angular) ray

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

 dt               = ( dr/( rho * ( rho + dr ) ) ) * t_prev * ( aesvt(j,1,ji_ray,jk_ray)/aesvt(j,2,ji_ray,jk_ray) )
 t                = dt + t_prev

 CALL eqstt_y( 1, j, ji_ray, jk_ray, rho + 0.5d0 * dr, 0.5d0 * ( t + t_prev ), ye, dum1, dum2, dpdt, dum3 )
 CALL eqstt_y( 2, j, ji_ray, jk_ray, rho + 0.5d0 * dr, 0.5d0 * ( t + t_prev ), ye, dum1, dum2, dedt, dum3)
 dt               = ( dr/( rho * ( rho + dr ) ) ) * 0.5d0 * ( t + t_prev ) * ( dpdt/dedt )
 t                = dt + t_prev

 CALL eqstt_y( 1, j, ji_ray, jk_ray,rho + 0.5d0 * dr, 0.5d0 * ( t + t_prev ), ye, dum1, dum2, dpdt, dum3 )
 CALL eqstt_y( 2, j, ji_ray, jk_ray,rho + 0.5d0 * dr, 0.5d0 * ( t + t_prev ), ye, dum1, dum2, dedt, dum3 )
 dt               = ( dr/( rho * ( rho + dr ) ) ) * 0.5d0 * ( t + t_prev ) * ( dpdt/dedt )
 t                = dt + t_prev

!........Done...........................................................

RETURN
END SUBROUTINE tgvndeye_dgen_y

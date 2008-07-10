SUBROUTINE w_cal( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, wgr, nx )
!-----------------------------------------------------------------------
!
!    File:         w_cal
!    Module:       w_cal
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/6/99
!
!    Purpose:
!      To compute the GR relativistic enthalpy.
!
!    Subprograms called:
!  eqstt_x     : interpolates quantities in the local EOS table
!
!    Input arguments:
!  jr_min      : inner zone for which calculation w is to be made
!  jr_max      : outer zone for which calculation w is to be made
!  ij_ray      : j-index of a radial ray
!  ik_ray      : k-index of a radial ray
!  rho(j)      : matter density of zone j (g/cm**3)
!  t(j)        : matter temperature of zone j (K)
!  ye(j)       : electron fraction of zone j
!
!    Output arguments:
!  wgr(j)      : relativistic enthalpy of zone j
!
!    Input arguments (common):
!  pq_x(j,ij_ray,ik_ray) : pseudoviscous pressure of zone j
!
!    Output arguments (common):
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, prb_cntl_module, shock_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : one, csqinv

USE edit_module, ONLY : nprint
USE prb_cntl_module, ONLY : irelhy
USE shock_module, ONLY : pq_x

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rho      ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: t        ! temperature (K)
REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: ye       ! electron fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx)   :: wgr      ! relativistic enthalpy

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index

REAL(KIND=double)                :: pmat          ! pressure (dynes cm^{-2})
REAL(KIND=double)                :: dpdd          ! d(pressure)/d(rho)
REAL(KIND=double)                :: dpdt          ! d(pressure)/d(temperature)
REAL(KIND=double)                :: dpdy          ! d(pressure)/d(electron fraction)
REAL(KIND=double)                :: umat          ! internal energy (ergs cm^{-3})
REAL(KIND=double)                :: dudd          ! d(internal energy)/d(rho)
REAL(KIND=double)                :: dudt          ! d(internal energy)/d(temperature)
REAL(KIND=double)                :: dudy          ! d(internal energy)/d(electron fraction)

  201 FORMAT (' irelhy =',i5,' not supported by w_cal')

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

GR: SELECT CASE (irelhy)

!-----------------------------------------------------------------------
!
!                 \\\\\ NONRELATIVISTIC HYDRO /////
!
!-----------------------------------------------------------------------

 CASE(0,1) GR

  DO j = jr_min,jr_max
    wgr(j)       = one
  END DO
  RETURN

!-----------------------------------------------------------------------
!
!                  \\\\\ RELATIVISTIC HYDRO /////
!
!-----------------------------------------------------------------------

 CASE(2) GR

  DO j = jr_min,jr_max
    CALL eqstt_x( 1, j, ij_ray, ik_ray, rho(j), t(j), ye(j), pmat, dpdd, dpdt, dpdy )
    CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye(j), umat, dudd, dudt, dudy )
    wgr(j)       = one + csqinv * ( umat + ( pmat + pq_x(j,ij_ray, ik_ray) )/rho(j) )
  END DO

!-----------------------------------------------------------------------
!  Stop if irelhy /= 0 and irelhy /= 1
!-----------------------------------------------------------------------

 CASE DEFAULT GR

  WRITE (nprint, 201) irelhy
  STOP

END SELECT GR

RETURN
END SUBROUTINE w_cal

SUBROUTINE gammaz_x( jr_min, jr_max, rho, t, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         gammaz_x
!    Module:       gammaz_x
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute the three adiabatic exponents, gam1(j,i_ray), gam2(j,i_ray),
!       and gam3(j,i_ray), for all radial zones j between jr_min and jr_max.
!       Here
!
!                       d(lnP)
!          gam1(j)   = --------    (constant s and ye)
!                      d(lnrho)
!
!          gam2(j)     d(lnP)
!        ----------- = ------      (constant s and ye)
!        gam2(j) - 1   d(lnt)
!
!                       d(lnt)
!        gam3(j) - 1 = --------    (constant s and ye)
!                      d(lnrho)
!
!    Input arguments:
!  jr_min           : minimum radial zone index.
!  jr_max           : maximum radial zone index.
!  rho              : shifted matter density array (g/cm**3).
!  t                : shifted matter matter temperature array (K).
!  ij_ray           : j-index of a radial ray
!  ik_ray           : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  aesv(j,i,i_ray)  : interpolated equation of state quantity i for radial zone j.
!  aesvd(j,i,i_ray) : d(aesv(j,i))/d(rho).
!  aesvd(j,i,i_ray  : d(aesv(j,i))/d(t).
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module,
!  eos_snc_x_module, mdl_cnfg_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx
USE numerical_module, ONLY : one, epsilon

USE eos_snc_x_module, ONLY : aesv, aesvd, aesvt, gam1, gam2, gam3

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min          ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max          ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t   ! shifted matter matter temperature array (K)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j               ! radial zone index

REAL(KIND=double)                :: gammat          ! temporary variable

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  gam3(j,ij_ray,ik_ray) = aesvt(j,1,ij_ray,ik_ray)                                              &
&                       / ( aesvt(j,2,ij_ray,ik_ray) * rho(j) + epsilon ) + one
  gam1(j,ij_ray,ik_ray) = ( t(j) * aesvt(j,1,ij_ray,ik_ray)                                     &
&                       * ( gam3(j,ij_ray,ik_ray) - one ) + rho(j) * aesvd(j,1,ij_ray,ik_ray) ) &
&                       / ( aesv(j,1,ij_ray,ik_ray) + epsilon )
  gammat                = ( rho(j) * aesvd(j,1,ij_ray,ik_ray)                                   &
&                       / ( gam3(j,ij_ray,ik_ray) - one ) + t(j) * aesvt(j,1,ij_ray,ik_ray) )   &
&                       / ( aesv(j,1,ij_ray,ik_ray) + epsilon )
  gam2(j,ij_ray,ik_ray) = gammat/( gammat - one )
END DO

RETURN
END SUBROUTINE gammaz_x

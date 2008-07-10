SUBROUTINE pseudo_x( j, ij_ray, ik_ray, rho_i, u, u_jm1, rho, r, r_jm1 )
!-----------------------------------------------------------------------
!
!    File:         pseudo_x
!    Module:       pseudo_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/14/01
!
!    Purpose:
!      To calculate the pseudoviscous pressure.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  j           : coordinate array index
!  ij_ray      : index denoting the j-index of a specific radial ray
!  ik_ray      : index denoting the k-index of a specific radial ray
!  rho_i       : matter density of zone j at time n (g/cm**3)
!  u           : velocity of zone j (cm)
!  rho         : matter density of zone j at time n+1 (g/cm**3)
!  r           : radius of zone j (cm)
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  ipq         : 0, pseudoviscosity (pq_x(j,ij_ray,ik_ray)) computed normally
!              : 1, pq_x(j,ij_ray,ik_ray)=0 unless u(jp) gt 0. for some jp
!              : 2, pq_x(j,ij_ray,ik_ray)=0
!
!    Output arguments (common):
!  pq_x(j,ij_ray,ik_ray) : pseudoviscous pressure
!
!    Include files:
!  kind_module, numerical_module
!  shock_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half

USE shock_module, ONLY : ipq, pq_x, pqr_x, q0_x

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! coordinate array index
INTEGER, INTENT(in)              :: ij_ray        ! iindex denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! iindex denoting the k-index of a specific radial ray

REAL(KIND=double), INTENT(in)    :: rho_i         ! initial density (cm^{-3})
REAL(KIND=double), INTENT(in)    :: u             ! velocity of the outer edge of zone j (cm s^{-1})
REAL(KIND=double), INTENT(in)    :: u_jm1         ! velocity of the inner edge of zone j (cm s^{-1})
REAL(KIND=double), INTENT(in)    :: rho           ! final density (cm^{-3})
REAL(KIND=double), INTENT(in)    :: r             ! radius of zone j (cm)
REAL(KIND=double), INTENT(in)    :: r_jm1         ! radius of zone j-1 (cm)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: lu            ! pseudoviscosity flag

REAL(KIND=double)                :: rajmh         ! zone-centered radius
REAL(KIND=double)                :: divv          ! divergence of fluid flow
REAL(KIND=double)                :: divvhm        ! divergence of fluid flow

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                     \\\\\ PSEUDOVISCOSITY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ipq: 0, pseudoviscosity (pq_x(j,ij_ray,ik_ray)) computed normally
!       1, pq_x(j,ij_ray,ik_ray)=0 unless u(jp) > 0 for some jp
!       2, pq_x(j,ij_ray,ik_ray)=0
!-----------------------------------------------------------------------

pq_x(j,ij_ray,ik_ray)   = zero

IF ( ipq == 2 ) RETURN

IF ( ipq == 1 ) THEN
  lu                    = .false.
  IF ( u >= 0.d0 ) THEN
    lu                  = .true.
  END IF ! u(l) .ge. 0
  IF ( .not. lu ) RETURN
END IF ! ipq = 1

rajmh                   = half * ( r + r_jm1 )
divv                    = ( r * r * u - r_jm1 * r_jm1 * u_jm1 )/( rajmh * rajmh )

IF ( j /= 2 ) divvhm    = rajmh * ( u/r - u_jm1/r_jm1 )
IF ( j == 2 ) divvhm    = zero

!-----------------------------------------------------------------------
!  pq_x(j,ij_ray,ik_ray) : pseudoviscous pressure at j-1/2
!
!  pq_x(j,ij_ray,ik_ray)=0 unless :
!  divv < 0  and  divvhm < 0  and  rho > rho_i
!-----------------------------------------------------------------------

pqr_x(j,ij_ray,ik_ray)  = pq_x(j,ij_ray,ik_ray)
pq_x(j,ij_ray,ik_ray)   = zero

IF ( rho > rho_i  .and.  divvhm <= zero  .and.  divv <= zero ) THEN
  pq_x(j,ij_ray,ik_ray) = q0_x(j) * rho * 4.d0 * (-divv) * (-divvhm)
END IF ! rho > rho, divvhm < 0, divv < 0

RETURN
END SUBROUTINE pseudo_x

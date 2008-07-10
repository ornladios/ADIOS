SUBROUTINE pseudo_z( j, ki_ray, kj_ray, rho_i, u, u_jm1, rho )
!-----------------------------------------------------------------------
!
!    File:         pseudo_z
!    Module:       pseudo_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/10/05
!
!    Purpose:
!      To calculate the pseudoviscous pressure.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  j           : coordinate array index
!  ki_ray      : x (radial) index of a specific z (azimuthal) ray
!  kj_ray      : y (azimuthal) index of a specific z (azimuthal) ray
!  rho_i       : matter density of zone j at time n (g/cm**3)
!  u           : velocity of zone j (cm)
!  rho         : matter density of zone j at time n+1 (g/cm**3)
!  r           : radius of zone j (cm)
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  ipq         : 0, pseudoviscosity (pq_z(j,kj_ray,ki_ray)) computed normally
!              : 1, pq_z(j,kj_ray,ki_ray)=0 unless u(jp) gt 0. for some jp
!              : 2, pq_z(j,kj_ray,ki_ray)=0
!
!    Output arguments (common):
!  pq_z(j,kj_ray,ki_ray) : pseudoviscous pressure
!
!    Include files:
!      kind_module, numerical_module
!      shock_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half

USE shock_module, ONLY : ipq, pq_z, pqr_z, q0_z

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! coordinate array index
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (azimuthal) index of a specific z (azimuthal) ray

REAL(KIND=double), INTENT(in)    :: rho_i         ! initial density (cm^{-3})
REAL(KIND=double), INTENT(in)    :: u             ! velocity of zone j (cm s^{-1})
REAL(KIND=double), INTENT(in)    :: u_jm1         ! velocity of zone j-1 (cm s^{-1})
REAL(KIND=double), INTENT(in)    :: rho           ! final density (cm^{-3})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)                :: divv          ! divergence of fluid flow
REAL(KIND=double)                :: divvhm        ! divergence of fluid flow

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Pseudoviscosity................................................
!.......................................................................
!.......................................................................

!-----------------------------------------------------------------------
!        ipq:      0, pseudoviscosity (pq_z(j,kj_ray,ki_ray)) computed normally
!                  1, pq_z(j,kj_ray,ki_ray)=0 unless u(jp) > 0 for some jp
!                  2, pq_z(j,kj_ray,ki_ray)=0
!-----------------------------------------------------------------------

pq_z(j,kj_ray,ki_ray)   = zero

IF ( ipq == 2 ) RETURN

divv                    = u - u_jm1
divvhm                  = u - u_jm1

!-----------------------------------------------------------------------
!        pq_z(j,kj_ray,ki_ray) :      pseudoviscous pressure at j-1/2
!        pq_z(j,kj_ray,ki_ray)=0 unless :
!          divv < 0  and  divvhm < 0  and  rho > rho_i
!-----------------------------------------------------------------------
pqr_z(j,kj_ray,ki_ray)  = pq_z(j,kj_ray,ki_ray)
pq_z(j,kj_ray,ki_ray)   = zero

IF ( rho > rho_i  .and.  divv <= zero ) THEN
  pq_z(j,kj_ray,ki_ray) = q0_z(j) * rho * 4.d0 * (-divv) * (-divvhm)
END IF ! rho > rho, divv < 0

RETURN
END SUBROUTINE pseudo_z

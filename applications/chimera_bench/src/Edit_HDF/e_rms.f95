SUBROUTINE e_rms( jr_min, jr_max, i_ray, i_ray_dim, nx, nnu, e_rms_stat, &
& e_rms_trns )
!-----------------------------------------------------------------------
!
!    File:         e_rms
!    Module:       e_rms
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/23/05
!
!    Purpose:
!      To compute the static and transport rms neutrino energies
!       as a function of radius and neutrino flavor.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min          : inner radial zone number
!  jr_max          : outer radial zone number
!  i_ray           : index denoting a specific radial ray
!  i_ray_dim       : number radial rays assigned to a processor
!  nx              : x-array extent
!  nnu             : neutrino flavor array extent
!
!    Output arguments:
!  e_rms_stat      : rms static neutrino energy
!  e_rms_trns      : rms transport neutino energy
!
!    Include files:
!  kind_module, numerical_module
!  nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon

USE nu_dist_module, ONLY : unu, dunu, psi0, psi1
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: jr_min        ! minimum radial index
INTEGER, INTENT(in)               :: jr_max        ! maximum radial index
INTEGER, INTENT(in)               :: i_ray         ! index denoting a specific radial ray
INTEGER, INTENT(in)               :: i_ray_dim     ! number radial rays assigned to a processor
INTEGER, INTENT(in)               :: nx            ! x-array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,nnu,i_ray_dim) :: e_rms_stat  ! SUM psi0 * w5dw/SUM w3ww
REAL(KIND=double), INTENT(out), DIMENSION(nx,nnu,i_ray_dim) :: e_rms_trns  ! SUM psi1 * w5dw/SUM w3ww

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index

REAL(KIND=double)                 :: w3dw          ! unu(j)**3 * dunu(j)
REAL(KIND=double)                 :: w5dw          ! unu(j)**5 * dunu(j)

REAL(KIND=double)                 :: psi0j3        ! psi0(j,k,n) * unu(j)**3 * dunu(j)
REAL(KIND=double)                 :: psi0j5        ! psi0(j,k,n) * unu(j)**5 * dunu(j)
REAL(KIND=double)                 :: psi1j3        ! psi1(j,k,n) * unu(j)**3 * dunu(j)
REAL(KIND=double)                 :: psi1j5        ! psi1(j,k,n) * unu(j)**5 * dunu(j)

REAL(KIND=double)                 :: enuv_stat2    ! psi0j5/psi0j3
REAL(KIND=double)                 :: enuv_trns2    ! psi1j5/psi1j3

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO j = jr_min,jr_max
    e_rms_stat(j,n,i_ray) = zero
    e_rms_trns(j,n,i_ray) = zero
  END DO
END DO

!-----------------------------------------------------------------------
!  Return if no neutrinos
!-----------------------------------------------------------------------

IF (nnugpmx == 0) RETURN

!-----------------------------------------------------------------------
!
!                \\\\\ NEUTRINO RMS ENERGIES /////
!
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO j = jr_min,jr_max

    psi0j3          = zero
    psi0j5          = zero
    psi1j3          = zero
    psi1j5          = zero

    IF ( nnugp(n) == 0 ) CYCLE

    DO k = 1,nnugp(n)

      w3dw          = unu(j,k)**3 * dunu(j,k)
      w5dw          = unu(j,k)**5 * dunu(j,k)

      psi0j3        = psi0j3 + w3dw  * psi0(j,k,n)
      psi0j5        = psi0j5 + w5dw  * psi0(j,k,n)
      psi1j3        = psi1j3 + w3dw  * psi1(j,k,n)
      psi1j5        = psi1j5 + w5dw  * psi1(j,k,n)

    END DO ! k

    enuv_stat2      = psi0j5 /( psi0j3 + epsilon )
    enuv_trns2      = psi1j5 /( psi1j3 + epsilon )

    e_rms_stat(j-jr_min+1,n,i_ray) = DSQRT( DABS(enuv_stat2) + epsilon )
    e_rms_trns(j-jr_min+1,n,i_ray) = DSQRT( DABS(enuv_trns2) + epsilon )

  END DO ! j
END DO ! n

RETURN
END SUBROUTINE e_rms

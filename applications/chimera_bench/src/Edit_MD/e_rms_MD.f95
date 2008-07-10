SUBROUTINE e_rms_MD( imin, imax, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& nx, nez, nnu, psi0, psi1, unu, dunu, e_rms_stat, e_rms_trns )
!-----------------------------------------------------------------------
!
!    File:         e_rms_MD
!    Module:       e_rms_MD
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
!  imin            : inner unshifted radial zone number
!  imax            : outer unshifted radial zone number
!  ij_ray          : j-index of a radial ray
!  ik_ray          : k-index of a radial ray
!  ij_ray_dim      : number of y-zones on a processor before swapping
!  ik_ray_dim      : number of z-zones on a processor before swapping
!  nx              : x-array extent
!  nez             : neutrino energy array extent
!  nnu             : neutrino flavor array extent
!  psi0            : zero moment of the neutrino distribution
!  psi1            : first moment of the neutrino distribution
!  unu             : zone-centered neutrino energy
!  dunu            : neutrino energy zone width
!
!    Output arguments:
!  e_rms_stat      : rms static neutrino energy
!  e_rms_trns      : rms transport neutino energy
!
!    Include files:
!  kind_module, numerical_module
!  nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon

USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: imin          ! inner unshifted radial zone number
INTEGER, INTENT(in)               :: imax          ! outer unshifted radial zone number
INTEGER, INTENT(in)               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)               :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)               :: nx            ! x-array extent
INTEGER, INTENT(in)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0       ! zero moment of f_neutrino
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1       ! zero moment of f_neutrino
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim)     :: unu        ! zone-centered neutrino energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,ij_ray_dim,ik_ray_dim)     :: dunu       ! neutrino energy zone width

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)    :: e_rms_stat ! SUM psi0 * w5dw/SUM w3ww
REAL(KIND=double), INTENT(out), DIMENSION(nx,nnu,ij_ray_dim,ik_ray_dim)    :: e_rms_trns ! SUM psi1 * w5dw/SUM w3ww

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: i             ! radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index

REAL(KIND=double)                 :: w3dw          ! unu(i)**3 * dunu(i)
REAL(KIND=double)                 :: w5dw          ! unu(i)**5 * dunu(i)

REAL(KIND=double)                 :: psi0j3        ! psi0(i,k,n) * unu(i)**3 * dunu(i)
REAL(KIND=double)                 :: psi0j5        ! psi0(i,k,n) * unu(i)**5 * dunu(i)
REAL(KIND=double)                 :: psi1j3        ! psi1(i,k,n) * unu(i)**3 * dunu(i)
REAL(KIND=double)                 :: psi1j5        ! psi1(i,k,n) * unu(i)**5 * dunu(i)

REAL(KIND=double)                 :: enuv_stat2    ! psi0j5/psi0j3
REAL(KIND=double)                 :: enuv_trns2    ! psi1j5/psi1j3

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

e_rms_stat(:,:,ij_ray,ik_ray) = zero
e_rms_trns(:,:,ij_ray,ik_ray) = zero

!-----------------------------------------------------------------------
!  Return if no neutrinos
!-----------------------------------------------------------------------

IF ( nnugpmx == 0 ) RETURN

!-----------------------------------------------------------------------
!
!                \\\\\ NEUTRINO RMS ENERGIES /////
!
!-----------------------------------------------------------------------

DO n = 1,nnu
  DO i = imin,imax

    psi0j3          = zero
    psi0j5          = zero
    psi1j3          = zero
    psi1j5          = zero

    IF ( nnugp(n) == 0 ) CYCLE

    DO k = 1,nnugp(n)

      w3dw          = unu(i,k,ij_ray,ik_ray)**3 * dunu(i,k,ij_ray,ik_ray)
      w5dw          = unu(i,k,ij_ray,ik_ray)**5 * dunu(i,k,ij_ray,ik_ray)

      psi0j3        = psi0j3 + w3dw  * psi0(i,k,n,ij_ray,ik_ray)
      psi0j5        = psi0j5 + w5dw  * psi0(i,k,n,ij_ray,ik_ray)
      psi1j3        = psi1j3 + w3dw  * psi1(i,k,n,ij_ray,ik_ray)
      psi1j5        = psi1j5 + w5dw  * psi1(i,k,n,ij_ray,ik_ray)

    END DO ! k

    enuv_stat2      = psi0j5 /( psi0j3 + epsilon )
    enuv_trns2      = psi1j5 /( psi1j3 + epsilon )

    e_rms_stat(i-imin+1,n,ij_ray,ik_ray) = DSQRT( DABS(enuv_stat2) + epsilon )
    e_rms_trns(i-imin+1,n,ij_ray,ik_ray) = DSQRT( DABS(enuv_trns2) + epsilon )

  END DO ! i
END DO ! n

RETURN
END SUBROUTINE e_rms_MD

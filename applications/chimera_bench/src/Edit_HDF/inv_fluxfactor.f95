SUBROUTINE inv_fluxfactor( jr_min, jr_max, i_ray, i_ray_dim, nx, nnu, &
& inv_fluxfact )
!-----------------------------------------------------------------------
!
!    File:         inv_fluxfactor
!    Module:       inv_fluxfactor
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/23/05
!
!    Purpose:
!      To compute the inverse flux factors as a function of radius
!       and neutrino flavor.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!   jr_min          : inner radial zone number
!   jr_max          : outer radial zone number
!   i_ray           : index denoting a specific radial ray
!   i_ray_dim       : number radial rays assigned to a processor
!   nx              : x-array extent
!   nnu             : neutrino flavor array extent
!
!    Output arguments:
!   inv_fluxfact    : inverse flux factor
!
!    Include files:
!  kind_module, numerical_module
!  mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, epsilon

USE mdl_cnfg_module, ONLY : r
USE nu_dist_module, ONLY : unue, dunue, psi0, psi1, rjmh
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

REAL(KIND=double), INTENT(out), DIMENSION(nx,nnu,i_ray_dim) :: inv_fluxfact ! inverse flux factor

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: j             ! radial zone index
INTEGER                           :: k             ! neutrino energy index
INTEGER                           :: n             ! neutrino flavor index

REAL(KIND=double)                 :: psi0j         ! psi0 interpolated to the zone edge
REAL(KIND=double)                 :: e_density     ! the neutrino energy density
REAL(KIND=double)                 :: e_flux        ! the neutrino flux

!........Initialize

DO n = 1,nnu
  DO j = jr_min,jr_max
    inv_fluxfact(j-jr_min+1,n,i_ray) = zero
  END DO
END DO

!........Return if no neutrinos

IF (nnugpmx == 0) RETURN

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO INVERSE FLUX FACTORS /////
!
!-----------------------------------------------------------------------

DO n = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO j = jr_min,jr_max
    e_density         = zero
    e_flux            = zero
    DO k = 1,nnugp(n)
      psi0j           = ( ( rjmh(j+1)**2 - r(j)**2 ) * psi0(j,k,n) &
&                     + ( r(j)**2 - rjmh(j)**2 ) * psi0(j+1,k,n) ) &
&                     / ( rjmh(j+1)**2 - rjmh(j)**2 )
      e_density       = e_density + unue(j,k)**3 * dunue(j,k) * psi0j
      e_flux          = e_flux + unue(j,k)**3 * dunue(j,k) * psi1(j,k,n)
    END DO ! K
    inv_fluxfact(j-jr_min+1,n,i_ray) = e_density/( e_flux + epsilon )
  END DO ! j
END DO ! n

RETURN
END SUBROUTINE inv_fluxfactor

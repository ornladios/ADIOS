SUBROUTINE radhyd_to_nu_energy_flux( ij_ray_dim, ik_ray_dim, nx, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nu_energy_flux
!    Module:       radhyd_to_nu_energy_flux
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To transfer variables from radial_ray_module, call the subroutines
!       to compute the total neutrino energy density and energy flux, and
!       transfer the results back to radial_ray_module
!
!    Subprograms called:
!  nu_energy  : computes the total neutrino energy density
!  nu_flux    : computes the total neutrino energy flux
!
!    Input arguments:
!  ij_ray_dim   : number of y-zones on a processor before swapping with y
!  ik_ray_dim   : number of z-zones on a processor before swapping with z
!  nx         : x_array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpi

USE radial_ray_module, ONLY : imin, imax, psi0_c, psi1_e, e_nu_c, f_nu_e
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                     :: ij_ray_dim ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)                     :: ik_ray_dim ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)                     :: nx         ! x-array extent
INTEGER, INTENT(in)                     :: nez        ! neutrino energy array extent
INTEGER, INTENT(in)                     :: nnu        ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                  :: jr_min    ! minimum shifted radial zone index
INTEGER                                  :: jr_max    ! maximum shifted radial zone index
INTEGER                                  :: ij_ray    ! j-index of a radial ray
INTEGER                                  :: ik_ray    ! k-index of a radial ray

REAL(KIND=double), DIMENSION(nx,nez,nnu) :: psi0      ! zero moment of the neutrino occupation probability
REAL(KIND=double), DIMENSION(nx,nez,nnu) :: psi1      ! first moment of the neutrino occupation probability
REAL(KIND=double), DIMENSION(nx)         :: u_nu      ! neutrino energy density (ergs cm^{-3})
REAL(KIND=double), DIMENSION(nx)         :: f_nu      ! neutrino energy flux (ergs cm^{-2} s^{-1})

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize shifted radial index boundaries
!-----------------------------------------------------------------------

jr_min                      = imin + 1
jr_max                      = imax + 1

!-----------------------------------------------------------------------
!
!              ||||| THE LOOP OVER THE RADIAL RAYS |||||
!              ||||| ON EACH PROCESSOR BEGINS HERE |||||
!
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Transfer variables from radial_ray_module to subroutine lists
!-----------------------------------------------------------------------

    psi0(jr_min:jr_max,:,:)   = psi0_c(imin:imax  ,:,:,ij_ray,ik_ray)
    psi1(imin:imax+1,:,:)     = psi1_e(imin:imax+1,:,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Compute neutrino energy densities and energy fluxes
!-----------------------------------------------------------------------

    CALL nu_energy( jr_min, jr_max, nx, nez, nnu, psi0, u_nu )
    CALL nu_flux( jr_min, jr_max, nx, nez, nnu, psi1, f_nu )

!-----------------------------------------------------------------------
!  Transfer neutrino energy densities and energy fluxes to
!   radial_ray_module
!-----------------------------------------------------------------------

    e_nu_c(imin:imax,ij_ray,ik_ray)   = u_nu(jr_min:jr_max)
    f_nu_e(imin:imax+1,ij_ray,ik_ray) = f_nu(imin:imax+1)

!-----------------------------------------------------------------------
!
!              ||||| THE LOOP OVER THE RADIAL RAYS |||||
!              |||||  ON EACH PROCESSOR ENDS HERE  |||||
!
!-----------------------------------------------------------------------

  END DO ! ij_ray
END DO ! ik_ray

RETURN
END SUBROUTINE radhyd_to_nu_energy_flux

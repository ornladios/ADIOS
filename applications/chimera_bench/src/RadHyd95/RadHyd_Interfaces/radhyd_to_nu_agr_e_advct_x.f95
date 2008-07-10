SUBROUTINE radhyd_to_nu_agr_e_advct_x( nx, ij_ray_dim, ik_ray_dim, nez, &
& nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nu_agr_e_advct_x
!    Module:       radhyd_to_nu_agr_e_advct_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the neutrino energy advection modules via subroutine
!       nu_energy_agr_advct_inout_x.
!
!    Subprograms called:
!  nu_energy_advct_inout_x : transfers variables to and from neutrino x-advection modules
!
!    Input arguments:
!  nx         : x-array extent
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  prb_cntl_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE prb_cntl_module, ONLY : ireltrns
USE radial_ray_module, ONLY : imin, imax, psi0_c, agr_c, agr_c_r

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ij_ray        ! j-index of a radial ray
INTEGER                          :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), DIMENSION(nx)         :: agrh  ! unshifted initial zone-centered lapse function
REAL(KIND=double), DIMENSION(nx)         :: agrah ! unshifted final zone-centered lapse function
REAL(KIND=double), DIMENSION(nx,nez,nnu) :: psi0  ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!  Return if ireltrns = 0
!-----------------------------------------------------------------------

IF ( ireltrns == 0 ) RETURN

!-----------------------------------------------------------------------
!            ||||| BEGIN LOOP OVER THE RADIAL RAYS |||||
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim

!-----------------------------------------------------------------------
!  Transfer neutrino e_advection lapse functions
!-----------------------------------------------------------------------

    agrh (imin:imax)        = agr_c_r(imin:imax,ij_ray,ik_ray)
    agrah(imin:imax)        = agr_c  (imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer neutrino distribution zero moment
!-----------------------------------------------------------------------

    psi0(imin:imax,:,:)     = psi0_c(imin:imax,:,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Do the energy advection
!-----------------------------------------------------------------------

    CALL nu_energy_agr_advct_inout_x( imin, imax, nx, nez, nnu, agrh,   &
&    agrah, psi0 )

!-----------------------------------------------------------------------
!  Transfer updated neutrino distribution zero moment back to
!   radial_ray_module
!-----------------------------------------------------------------------

    psi0_c(imin:imax,:,:,ij_ray,ik_ray)                                 &
&                           = psi0(imin:imax,:,:)

!-----------------------------------------------------------------------
!             ||||| END LOOP OVER THE RADIAL RAYS |||||
!-----------------------------------------------------------------------

  END DO ! ij_ray
END DO ! ik_ray

RETURN
END SUBROUTINE radhyd_to_nu_agr_e_advct_x

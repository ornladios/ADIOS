SUBROUTINE radhyd_to_transport( nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, &
& nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_transport
!    Module:       radhyd_to_transport
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the x-array neutrino transport modules via subroutine
!       nu_transport_inout.
!
!    Input arguments:
!  nx                 : x-array extent
!  ij_ray_dim         : number of y-zones on a processor before swapping with y
!  ik_ray_dim         : number of z-zones on a processor before swapping with z
!  ij_ray             : index denoting the j-index of a specific radial ray
!  ik_ray             : index denoting the k-index of a specific radial ray
!  nez                : neutrino energy array extent
!  nnu                : neutrino flavor array extent
!  nnc                : composition array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  nu_transport_inout : directs the computation of the neutrino transport
!
!    Include files:
!  kind_module
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : imin, imax, nprint, ncycle, dtime=>dtnph, &
& rho_c, t_c, ye_c, rhobar, x_ei, u_c, agr_e, agr_c, psi0_c, psi1_e,    &
& rhs1_c, dc_e, dt, jdt, dtnph_trans, dtime_trans, xn_c

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc           ! composition array extent

!-----------------------------------------------------------------------
!  Transfer variables to and from mgfld
!-----------------------------------------------------------------------

CALL nu_transport_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim,     &
& ik_ray_dim, nez, nnu, nnc, nprint, rho_c, t_c, ye_c, rhobar, x_ei,     &
& u_c, agr_e, agr_c, psi0_c, psi1_e, rhs1_c, dc_e, dt, jdt, dtnph_trans, &
& dtime_trans(ij_ray,ik_ray), xn_c )

RETURN
END SUBROUTINE radhyd_to_transport

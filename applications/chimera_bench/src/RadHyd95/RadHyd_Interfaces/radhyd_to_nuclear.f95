SUBROUTINE radhyd_to_nuclear( nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, &
& nnc ) 
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nuclear
!    Module:       radhyd_to_nuclear
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the nuclear burn modules via subroutine nuc_burn_inout.
!
!    Input arguments:
!  nx             : x-array extent
!  ij_ray_dim     : number of y-zones on a processor before swapping
!  ik_ray_dim     : number of z-zones on a processor before swapping
!  ij_ray         : index denoting the j-index of a specific radial ray
!  ik_ray         : index denoting the k-index of a specific radial ray
!  nnc            : composition array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  nuc_burn_inout : directs the computation of the nuclear evolution
!
!    Include files:
!  kind_module
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : imin, imax, nprint, ncycle, dtime=>dtnph, &
& time, rho_c, t_c, ye_c, xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, &
& nse_c, dgrid, tgrid, ygrid, rhoes, idty, dt, jdt

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: ij_ray        ! iindex denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! iindex denoting the k-index of a specific radial ray
INTEGER, INTENT(in)              :: nnc           ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls            ! inner logical composition zone
INTEGER                          :: le            ! outer logical composition zone

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                    = 1
le                    = nnc

!-----------------------------------------------------------------------
!  Transfer variables to nuclear arrays
!-----------------------------------------------------------------------

CALL nuc_burn_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& ls, le, nnc, dgrid, tgrid, ygrid, rhoes, idty, nprint, ncycle, dtime, &
& time, rho_c, t_c, ye_c, xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, &
& nse_c, dt, jdt )

RETURN
END SUBROUTINE radhyd_to_nuclear

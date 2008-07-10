SUBROUTINE nse_set( nx, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         nse_set
!    Module:       nse_set
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/31/05
!
!    Purpose:
!      To set nse, the nuclear statistical equilibrium flag, for all
!       x-array zones beyond the outer edge of the model, to the value
!       of nse at the outer edge of the model.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x-array array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  eos_snc_x_module, nucbrn_module, radial_ray_module, 
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE eos_snc_x_module, ONLY : nse_e=>nse
USE nucbrn_module, ONLY : nse
USE radial_ray_module, ONLY : imax, nse_c

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx           ! energy array dimension
INTEGER, INTENT(in)              :: ij_ray_dim   ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim   ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ij_ray       ! j-index of a radial ray
INTEGER                          :: ik_ray       ! k-index of a radial ray

!-----------------------------------------------------------------------
!
!              \\\\\ SET NES FALG FOR OUTER GHOSTS /////
!
!-----------------------------------------------------------------------

DO ik_ray = 1,ik_ray_dim
  DO ij_ray = 1,ij_ray_dim
    nse_c(imax+1:nx,ij_ray,ik_ray) = nse_c(imax,ij_ray,ik_ray)
    nse_e(imax+1:nx,ij_ray,ik_ray) = nse_c(imax,ij_ray,ik_ray)
  END DO ! ij_ray = 1,ij_ray_dim
END DO ! ik_ray = 1,ik_ray_dim

nse(imax+1:nx)                     = nse_c(imax,1,1)

RETURN
END SUBROUTINE nse_set

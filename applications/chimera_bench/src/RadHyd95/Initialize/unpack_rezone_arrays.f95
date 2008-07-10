SUBROUTINE unpack_rezone_arrays( l_rezone_data, d_rezone_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_rezone_arrays
!    Module:       unpack_rezone_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/27/04
!
!    Purpose:
!      To unpack the rezone array data defining the y anz z coordinates
!       of the problem.
!
!    Subprograms called:
!
!    Input arguments:
!  l_rezone_data : integer array of reaone keys
!  d_rezone_data : 64 bit real array of reaone keys
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module
!  evh1_global, evh1_sweep, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : ny, nz

USE evh1_global, ONLY : lagrangian, courant
USE evh1_sweep, ONLY : xmin_s=>xmin, xmax_s=>xmax, ymin_s=>ymin, ymax_s=>ymax, &
& zmin_s=>ymin, zmax_s=>ymax
USE radial_ray_module, ONLY : y_ei, y_ci, dy_ci, z_ei, z_ci, dz_ci, xmin, &
& xmax, ymin, ymax, zmin, zmax, xmin_i, xmax_i, ymin_i, ymax_i, zmin_i, zmax_i

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

LOGICAL, INTENT(inout), DIMENSION(1)           :: l_rezone_data ! integer array of reaone keys

REAL(KIND=double), INTENT(inout), DIMENSION(20+3*ny+3*nz+2) :: d_rezone_data ! 64 bit real array of rezone data

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                                        :: i          ! do inddex

!-----------------------------------------------------------------------
!
!                  \\\\\ UNPACK REZONE KEYS /////
!
!-----------------------------------------------------------------------

lagrangian          = l_rezone_data(1)

xmin                = d_rezone_data(1)
xmax                = d_rezone_data(2)
ymin                = d_rezone_data(3)
ymax                = d_rezone_data(4)
zmin                = d_rezone_data(5)
zmax                = d_rezone_data(6)
xmin_s              = d_rezone_data(1)
xmax_s              = d_rezone_data(2)
ymin_s              = d_rezone_data(3)
ymax_s              = d_rezone_data(4)
zmin_s              = d_rezone_data(5)
zmax_s              = d_rezone_data(6)
xmin_i              = d_rezone_data(7)
xmax_i              = d_rezone_data(8)
ymin_i              = d_rezone_data(9)
ymax_i              = d_rezone_data(10)
zmin_i              = d_rezone_data(11)
zmax_i              = d_rezone_data(12)
courant             = d_rezone_data(13)

DO i = 1,ny
  y_ei(i)           = d_rezone_data(20+i)
  y_ci(i)           = d_rezone_data(20+ny+1+i)
  dy_ci(i)          = d_rezone_data(20+2*ny+1+i)
END DO
y_ei(ny+1)          = d_rezone_data(20+ny+1)

DO i = 1,nz
  z_ei(i)           = d_rezone_data(20+3*ny+1+i)
  z_ci(i)           = d_rezone_data(20+3*ny+nz+2+i)
  dz_ci(i)          = d_rezone_data(20+3*ny+2*nz+2+i)
END DO
z_ei(nz+1)          = d_rezone_data(20+3*ny+nz+2)

RETURN
END SUBROUTINE unpack_rezone_arrays

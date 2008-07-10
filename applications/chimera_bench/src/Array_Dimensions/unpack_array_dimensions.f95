SUBROUTINE unpack_array_dimenisons( c_path_data, n_dim_data, nx, ny, nz, &
& nez, nnu, nnc, n_proc, n_proc_y, n_proc_z, ij_ray_dim, ik_ray_dim, &
& j_ray_dim, k_ray_dim, data_path, log_path, reset_path )
!-----------------------------------------------------------------------
!
!    File:         unpack_array_dimenisons
!    Module:       unpack_array_dimenisons
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/02/04
!
!    Purpose:
!      To unpack the array dimensions defining the extents of the arrays.
!
!    Subprograms called:
!
!    Input arguments:
!  c_path_data : path to the output data directories
!  n_dim_data  : array dimension data
!
!    Output arguments:
!  nx          : x-array dimension
!  ny          : y (angular array) dimension
!  nz          : z dimension
!  nez         : neutrino energy array dimension
!  nnu         : neutrino flavor array dimension
!  nnc         : composition array dimension
!  n_proc      : number of processors assigned to the run
!  n_proc_y    : number of processors assigned to the y-zones
!  n_proc_z    : number of processors assigned to the z-zones
!  ij_ray_dim  : number of y-zones on a processor before swapping
!  ik_ray_dim  : number of z-zones on a processor before swapping
!  j_ray_dim   : number of radial zones on a processor after swapping with y
!  k_ray_dim   : number of radial zones on a processor after swapping with z
!  data_path   : path to the output data directories
!  log_path    : path to the simulation log (blank writes to the screen)
!  reset_path  : path to write the restart keys file, reset.d
!
!    Include files:
!        none
!
!-----------------------------------------------------------------------

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER(len=128), INTENT(in), DIMENSION(3) :: c_path_data ! character array of data path

INTEGER, INTENT(in), DIMENSION(20)           :: n_dim_data  ! array dimension data

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER (len = 128), INTENT(out) :: data_path     ! path to the output data directories
CHARACTER (len = 128), INTENT(out) :: log_path      ! path to the simulation log
CHARACTER (len = 128), INTENT(out) :: reset_path    ! path to write the restart keys file, reset.d

INTEGER, INTENT(out)               :: nx            ! x-array extent
INTEGER, INTENT(out)               :: ny            ! y-array extent
INTEGER, INTENT(out)               :: nz            ! z-array extent
INTEGER, INTENT(out)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(out)               :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(out)               :: nnc           ! composition array extent
INTEGER, INTENT(out)               :: n_proc        ! number of processors assigned to the run
INTEGER, INTENT(out)               :: n_proc_y      ! number of processors assigned to the y-zones
INTEGER, INTENT(out)               :: n_proc_z      ! number of processors assigned to the z-zones
INTEGER, INTENT(out)               :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(out)               :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(out)               :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(out)               :: k_ray_dim     ! number of radial zones on a processor after swapping with z
 
!-----------------------------------------------------------------------
!
!               \\\\\ UNPACK ARRAY DIMENSIONS /////
!
!-----------------------------------------------------------------------

nx             = n_dim_data(1)
ny             = n_dim_data(2)
nz             = n_dim_data(3)
nez            = n_dim_data(4)
nnu            = n_dim_data(5)
nnc            = n_dim_data(6)
n_proc         = n_dim_data(7)
n_proc_y       = n_dim_data(8)
n_proc_z       = n_dim_data(9)
ij_ray_dim     = n_dim_data(10)
ik_ray_dim     = n_dim_data(11)
j_ray_dim      = n_dim_data(12)
k_ray_dim      = n_dim_data(13)

data_path      = c_path_data(1)
log_path       = c_path_data(2)
reset_path     = c_path_data(3)

RETURN
END SUBROUTINE unpack_array_dimenisons

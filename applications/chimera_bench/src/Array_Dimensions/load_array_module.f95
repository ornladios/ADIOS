SUBROUTINE load_array_module( nxp, nyp, nzp, nezp, nnup, nncp, n_procp, &
& n_proc_yp, n_proc_zp, n_xray_yp, n_xray_zp, n_yrayp, n_zrazp, max_12p )
!-----------------------------------------------------------------------
!
!    File:         load_array_module
!    Module:       load_array_module
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/02/04
!
!    Purpose:
!        To load array dimensions into array_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nxp       : x-array extent
!  nyp       : y-array extent
!  nzp       : z-array extent
!  nezp      : neutrino energy array extent
!  nnup      : neutrino flavor array extent
!  nncp      : composition array extent
!  n_procp   : number of processors assigned to a run
!  n_proc_yp : number of processors assigned to the y-zones
!  n_proc_zp : number of processors assigned to the z-zones
!  n_xray_yp : number of y-zones on a processor before swapping with y
!  n_xray_yp : number of z-zones on a processor before swapping with z
!  n_yrayp   : number of radial zones on a processor after swapping with y
!  n_zrayp   : number of radial zones on a processor after swapping with z
!
!    Output arguments:
! max_12p    : max(nx,ny,nz)+12
!
!    Include files:
!  array_module
!
!-----------------------------------------------------------------------

USE array_module, ONLY : nx, ny, nz, nez, nnu, nnc, nezp1, nez2, nez2p, &
& max_12, n_proc, n_proc_y, n_proc_z, ij_ray_dim, ik_ray_dim, j_ray_dim, &
& k_ray_dim

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)    :: nxp           ! x-array extent
INTEGER, INTENT(in)    :: nyp           ! y-array extent
INTEGER, INTENT(in)    :: nzp           ! z-array extent
INTEGER, INTENT(in)    :: nezp          ! neutrino energy array extent
INTEGER, INTENT(in)    :: nnup          ! neutrino flavor array extent
INTEGER, INTENT(in)    :: nncp          ! composition array extent
INTEGER, INTENT(in)    :: n_procp       ! number of processors assigned to a run
INTEGER, INTENT(in)    :: n_proc_yp     ! number of processors assigned to the y-zones
INTEGER, INTENT(in)    :: n_proc_zp     ! number of processors assigned to the z-zones
INTEGER, INTENT(in)    :: n_xray_yp     ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)    :: n_xray_zp     ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)    :: n_yrayp       ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)    :: n_zrazp       ! number of radial zones on a processor after swapping with z

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out)   :: max_12p       ! max(nx,ny,nz)+12

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                   \\\\\ ARRAY DIMENSIONS /////
!
!-----------------------------------------------------------------------

nx             = nxp
ny             = nyp
nz             = nzp
nez            = nezp
nnu            = nnup
nnc            = nncp
nezp1          = nez + 1
nez2           = nez + 2
nez2p          = 2 * nez2 + 1
max_12         = MAX(nx, ny, nz) + 12
max_12p        = max_12
n_proc         = n_procp
n_proc_y       = n_proc_yp
n_proc_z       = n_proc_zp
ij_ray_dim     = n_xray_yp
ik_ray_dim     = n_xray_zp
j_ray_dim      = n_yrayp
k_ray_dim      = n_zrazp

RETURN
END SUBROUTINE load_array_module

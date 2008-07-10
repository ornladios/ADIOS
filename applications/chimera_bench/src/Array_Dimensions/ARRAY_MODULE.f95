!-----------------------------------------------------------------------
!    Module:       array_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!
!    Specification of the array dimensions
!-----------------------------------------------------------------------

MODULE array_module
SAVE

!-----------------------------------------------------------------------
!  nx, ny, nz, nez, nnu, nnc
!-----------------------------------------------------------------------
!  nx  : is the maximum number of radial zones
!  ny  : is the maximum number of y-zones
!  nz  : is the maximum number of z-zones
!  nez : is the maximum number of energy zones
!  nnu : is the maximum number of neutrino degrees of freedom
!  nnc : the maximum number of nuclei for matter not in nuclear statistical
!   equilibrium
!-----------------------------------------------------------------------

INTEGER                    :: nx
INTEGER                    :: ny
INTEGER                    :: nz
INTEGER                    :: nez
INTEGER                    :: nnu
INTEGER                    :: nnc

!-----------------------------------------------------------------------
!  nezp1, nez2, nez2p
!-----------------------------------------------------------------------

INTEGER                    :: nezp1
INTEGER                    :: nez2
INTEGER                    :: nez2p
INTEGER                    :: max_12

!-----------------------------------------------------------------------
!  n_proc
!-----------------------------------------------------------------------
!  n_proc     : the number of processors assigned to the run
!  n_proc_y   : the number of processors assigned to the y-zones
!  n_proc_z   : the number of processors assigned to the z-zones
!  ij_ray_dim : the number of y-zones on a processor before swapping with y
!  ik_ray_dim : the number of z-zones on a processor before swapping with z
!  j_ray_dim  : the number of radial zones on a processor after swapping
!                with y
!  k_ray_dim  : the number of radial zones on a processor after swapping
!                with z
!-----------------------------------------------------------------------

INTEGER                    :: n_proc
INTEGER                    :: n_proc_y
INTEGER                    :: n_proc_z
INTEGER                    :: ij_ray_dim
INTEGER                    :: ik_ray_dim
INTEGER                    :: j_ray_dim
INTEGER                    :: k_ray_dim

END module array_module

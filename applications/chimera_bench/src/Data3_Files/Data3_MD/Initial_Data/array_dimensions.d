

!-----------------------------------------------------------------------
!         RadHyd array dimensions
!
!      nx     : x-array extent
!      ny     : y-array extent
!      nz     : z-array extent
!      nez    : neutrino energy array extent
!      nnu    : neutrino flavor array extent
!      nnc    : abundance array extent
!      n_proc : number of processors assigned to the run
!
!-----------------------------------------------------------------------

nx                         300                                          nx
ny                          16                                          ny
nz                           4                                          nz
nez                         20                                          nez
nnu                          4                                          nnu
nnc                         17                                          nnc
n_proc                      16                                          n_proc
n_proc                       8                                          n_proc_y
n_proc                       2                                          n_proc_z

!-----------------------------------------------------------------------
!      data_path: the path directing output to the data directories.
!-----------------------------------------------------------------------

data_p    /scratch/scratchdirs/bruenn/Data3_3D

!-----------------------------------------------------------------------
!      log_path : the path directing output to the simulation log.
!-----------------------------------------------------------------------

log_p     /scratch/scratchdirs/bruenn/Data3_3D


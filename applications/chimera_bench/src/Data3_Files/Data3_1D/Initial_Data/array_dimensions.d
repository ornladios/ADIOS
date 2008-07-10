

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
ny                           1                                          ny
nz                           1                                          nz
nez                         20                                          nez
nnu                          4                                          nnu
nnc                         17                                          nnc
n_proc                       1                                          n_proc
n_proc                       1                                          n_proc_y
n_proc                       1                                          n_proc_z

!-----------------------------------------------------------------------
!      data_path: the path directing output to the data directories.
!-----------------------------------------------------------------------

data_p    /scratch/scratchdirs/bruenn/Data3_3D_1

!-----------------------------------------------------------------------
!      log_path : the path directing output to the simulation log.
!-----------------------------------------------------------------------

log_p     /scratch/scratchdirs/bruenn/Data3_3D_1


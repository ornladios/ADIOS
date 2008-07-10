

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
ny                           8                                          ny
nz                           8                                          nz
nez                         20                                          nez
nnu                          4                                          nnu
nnc                         17                                          nnc
n_proc                      64                                          n_proc
n_proc                       8                                          n_proc_y
n_proc                       8                                          n_proc_z

!-----------------------------------------------------------------------
!      data_path: the path directing output to the data directories.
!-----------------------------------------------------------------------

data_p    /tmp/work/zf2/RadHyd3D_adios/Execute/build/Data3

!-----------------------------------------------------------------------
!      log_path : the path directing output to the simulation log.
!-----------------------------------------------------------------------

log_p     /tmp/work/zf2/RadHyd3D_adios/Execute/build/Data3 

!-----------------------------------------------------------------------
!      reset_path : the path and directory in which to write the
!       restart keys file 'reset.d'
!      Set to Data3/Initial_Data to restart from the latest temp
!       restart file.
!-----------------------------------------------------------------------

rset_p    Data3/Initial_Data


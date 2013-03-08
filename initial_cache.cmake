set(BUILD_WRITE ON CACHE BOOL "")

set(BUILD_FORTRAN OFF CACHE BOOL "")

#set(MXML_DIR "/sw/analysis-x64/mxml/2.6/centos5.5_pgi10.9" CACHE FILEPATH "path to mxml dir")

set(DATATAP OFF CACHE BOOL "")

set(DATASPACES OFF CACHE BOOL "")

set(DIMES OFF CACHE BOOL "")

set(NSSI OFF CACHE BOOL "")

set(NC4PAR OFF CACHE BOOL "")

set(PHDF5 ON CACHE BOOL "")
# phdf5: /sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4
set(PHDF5_DIR "/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4" CACHE FILEPATH "path to parallel hdf5 dir")
set(PHDF5_FLAGS "-I/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4/include" CACHE FILEPATH "flags to use parallel hdf5")
set(PHDF5_LIBS "-L/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4/lib -lhdf5_hl -lhdf5 -L/sw/sith/szip/2.1/centos5.5_pgi10.9/lib -lsz -lz -lm" CACHE STRING "parallel hdf5")

set(PORTALS OFF CACHE BOOL "")

set(INFINIBAND OFF CACHE BOOL "")

set(CRAY_PMI OFF CACHE BOOL "")

set(CRAY_UGNI OFF CACHE BOOL "")

set(NETCDF OFF CACHE BOOL "")

set(HDF5 ON CACHE BOOL "")
set(HDF5_DIR "/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9" CACHE FILEPATH "path to hdf5 dir")
set(HDF5_FLAGS "-I/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9/include" CACHE FILEPATH "flags to use suquential hdf5")
set(HDF5_LIBS "-L/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9/lib -lhdf5_hl -lhdf5 -L/sw/sith/szip/2.1/centos5.5_pgi10.9/lib -lsz -lz -lm" CACHE FILEPATH "sequential hdf5")

set(DMALLOC OFF CACHE BOOL "")

set(LUSTRE OFF CACHE BOOL "")
set(LUSTRE_DIR "" CACHE FILEPATH "path to lustre dir")


set(BUILD_WRITE ON CACHE BOOL "")
set(BUILD_FORTRAN OFF CACHE BOOL "")

set(MXML_DIR "" CACHE FILEPATH "path to mxml dir")

set(PHDF5 ON CACHE BOOL "")
# phdf5: /sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4
set(PHDF5_FLAGS "-I/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4/include" CACHE FILEPATH "flags to use parallel hdf5")
set(PHDF5_LIBS "-L/sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4/lib -lhdf5_hl -lhdf5 -L/sw/sith/szip/2.1/centos5.5_pgi10.9/lib -lsz -lz -lm" CACHE STRING "parallel hdf5"

set(HDF5 OFF CACHE BOOL "")
set(HDF5_DIR "" CACHE FILEPATH "path to hdf5 dir")

set(LUSTRE OFF CACHE BOOL "")
set(LUSTRE_DIR "" CACHE FILEPATH "path to lustre dir")


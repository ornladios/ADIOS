if(DEFINED ENV{BUILD_WRITE})
  set(BUILD_WRITE $ENV{BUILD_WRITE})
  if(BUILD_WRITE)
    set(BUILD_WRITE ON CACHE BOOL "")
  else()
    set(BUILD_WRITE OFF CACHE BOOL "")
  endif()
else()
  set(BUILD_WRITE OFF CACHE BOOL "")
endif()

if(DEFINED ENV{BUILD_FORTRAN})
  set(BUILD_FORTRAN $ENV{BUILD_FORTRAN})
  if(BUILD_FORTRAN)
    set(BUILD_FORTRAN ON CACHE BOOL "")
  else()
    set(BUILD_FORTRAN OFF CACHE BOOL "")
  endif()
else()
  set(BUILD_FORTRAN OFF CACHE BOOL "")
endif()

set(MACRODEFFLAG $ENV{MACRODEFFLAG})

set(MXML_DIR $ENV{MXML_DIR}  CACHE FILEPATH "path to mxml dir")

set(DATATAP OFF CACHE BOOL "")

set(DIMES OFF CACHE BOOL "")

set(NSSI OFF CACHE BOOL "")

set(NC4PAR OFF CACHE BOOL "")

if(NC4PAR)
  set(NC4PAR_DIR "$ENV{PAR_NC_DIR}" CACHE FILEPATH "path to parallel NETCDF dir")
endif()

set(PHDF5 ON CACHE BOOL "")

if(PHDF5)
  set(PHDF5_DIR $ENV{PAR_HDF5_DIR} CACHE FILEPATH "path to parallel hdf5 dir")
  set(PHDF5_FLAGS "-I${PHDF5_DIR}/include" CACHE FILEPATH "flags to use parallel hdf5")
  set(PHDF5_LIBS $ENV{PAR_HDF5_CLIB} CACHE STRING "parallel hdf5")
endif()

set(PORTALS OFF CACHE BOOL "")

if(DEFINED ENV{DATASPACES_INCDIR})
  set(DATASPACES_INCDIR $ENV{DATASPACES_INCDIR} CACHE INTERNAL "Internal variable")
  if(DATASPACES_INCDIR)
    if(DEFINED $ENV{DATASPACES_LIBDIR})
      set(DATASPACES_LIBDIR $ENV{DATASPACES_LIBDIR} CACHE INTERNAL "Internal variable")
        if(DATASPACES_LIBDIR)
          set(DATASPACES ON CACHE BOOL "")
        endif()
    endif()
  endif()
elseif(DEFINED ENV{DATASPACES_DIR})
  set(DATASPACES ON CACHE BOOL "")
  set(DATASPACES_DIR $ENV{DATASPACES_DIR} CACHE FILEPATH "path to dataspaces dir")
#  set(DATASPACES_LIBS "-L$ENV{DATASPACES_DIR}/libs -ldspaces -ldscommon -ldar" CACHE STRING "dataspace libs")
elseif(DEFINED ENV{DATASPACES})
  set(DATASPACES ON CACHE BOOL "")
  set(DATASPACES_DIR $ENV{DATASPACES} CACHE FILEPATH "path to dataspaces dir")
else()
  set(DATASPACES OFF CACHE BOOL "")
endif()


set(CRAY_PMI OFF CACHE BOOL "")

set(CRAY_UGNI OFF CACHE BOOL "")

set(NETCDF OFF CACHE BOOL "")

if(NETCDF)
  set(NETCDF_DIR "$ENV{SEQ_NC_DIR}" CACHE FILEPATH "path to suquential NETCDF dir")
endif()

set(HDF5 ON CACHE BOOL "")

if(HDF5)
  set(HDF5_DIR "$ENV{SEQ_HDF5_DIR}" CACHE FILEPATH "path to hdf5 dir")
  set(HDF5_FLAGS "-I${HDF5_DIR}/include" CACHE FILEPATH "flags to use suquential hdf5")
  set(HDF5_LIBS "$ENV{SEQ_HDF5_CLIB}" CACHE FILEPATH "sequential hdf5")
endif()

set(DMALLOC OFF CACHE BOOL "")

set(LUSTRE OFF CACHE BOOL "")
set(LUSTRE_DIR "" CACHE FILEPATH "path to lustre dir")


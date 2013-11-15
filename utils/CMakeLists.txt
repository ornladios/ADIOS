add_subdirectory(gpp)
add_subdirectory(bpdump)
add_subdirectory(bp2ascii)
add_subdirectory(bpsplit)
add_subdirectory(bpls)
add_subdirectory(bpdiff)

if(BUILD_WRITE)
  add_subdirectory(adios_lint)
  add_subdirectory(bp2bp)
endif(BUILD_WRITE)

if(HAVE_HDF5)
  add_subdirectory(bp2h5)
endif(HAVE_HDF5)

if(HAVE_NETCDF)
  add_subdirectory(bp2ncd)
endif(HAVE_NETCDF)

add_subdirectory(skel)
add_subdirectory(skeldump)

#find_library(ADIOS_INTERNAL_NOMPI_LIB adios_internal_nompi ${PROJECT_BINARY_DIR}/src)
#if(${ADIOS_INTERNAL_NOMPI_LIB}-NOTFOUND)
#  message(STATUS "cannot find adios_internal_nompi.a in ${PROJECT_SOURCE_DIR}/utils")
#endif(${ADIOS_INTERNAL_NOMPI_LIB}-NOTFOUND)

#find_library(ADIOSREAD_NOMPI_LIB adiosread_nompi ${PROJECT_BINARY_DIR}/src)
#if(${ADIOSREAD_NOMPI_LIB}-NOTFOUND)
#  message(STATUS "cannot find adiosread_nompi.a in ${PROJECT_SOURCE_DIR}/utils")
#endif(${ADIOSREAD_NOMPI_LIB}-NOTFOUND)

#find_library(ADIOS_LIB adios ${PROJECT_BINARY_DIR}/src)
#if(${ADIOS_LIB}-NOTFOUND)
#  message(STATUS "cannot find adios.a in ${PROJECT_SOURCE_DIR}/utils makefile")
#endif(${ADIOS_LIB}-NOTFOUND)

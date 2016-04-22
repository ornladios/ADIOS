add_subdirectory(gpp)
add_subdirectory(bpdump)
add_subdirectory(bp2ascii)
add_subdirectory(bpsplit)
add_subdirectory(bpls)
add_subdirectory(bpdiff)
add_subdirectory(bpmeta)
add_subdirectory(bprecover)
add_subdirectory(adios_list_methods)

if(BUILD_WRITE)
  add_subdirectory(adios_lint)
  add_subdirectory(bp2bp)
  set_property(TARGET adios_lint PROPERTY LINKER_LANGUAGE CXX)
  set_property(TARGET bp2bp PROPERTY LINKER_LANGUAGE CXX)
endif(BUILD_WRITE)

if(HAVE_HDF5)
  add_subdirectory(bp2h5)
  set_property(TARGET bp2h5 PROPERTY LINKER_LANGUAGE CXX)
endif(HAVE_HDF5)

if(HAVE_NETCDF)
  add_subdirectory(bp2ncd)
  set_property(TARGET bp2ncd PROPERTY LINKER_LANGUAGE CXX)
  set_property(TARGET fastbit PROPERTY LINKER_LANGUAGE CXX)
endif(HAVE_NETCDF)

if(HAVE_FASTBIT)
  add_subdirectory(fastbit)
endif()

add_subdirectory(skel)
add_subdirectory(skeldump)

set_property(TARGET bprecover PROPERTY LINKER_LANGUAGE CXX)

set_property(TARGET skel PROPERTY LINKER_LANGUAGE CXX)
set_property(TARGET skeldump PROPERTY LINKER_LANGUAGE CXX)


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

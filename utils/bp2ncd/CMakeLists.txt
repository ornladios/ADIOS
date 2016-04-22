include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/utils/bp2ncd)
include_directories(${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/src  ${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/utils/bp2ncd)

add_executable(bp2ncd bp2ncd.c)

if(USE_PARALLEL_COMPILER)
  set(bp2ncd_CPPFLAGS "${NETCDF_CPPFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")
  set(bp2ncd_CFLAGS "${NETCDF_CFLAGS} ${ADIOSLIB_CFLAGS}")
  target_link_libraries(bp2ncd adios ${ADIOSLIB_LDADD})
  set_target_properties(bp2ncd PROPERTIES COMPILE_FLAGS "${bp2ncd_CPPFLAGS} ${bp2ncd_CFLAGS}")
else(USE_PARALLEL_COMPILER)
  set(bp2ncd_CPPFLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")
  target_link_libraries(bp2ncd adios_internal_nompi ${ADIOSLIB_INT_LDADD})
  set_target_properties(bp2ncd PROPERTIES COMPILE_FLAGS "${bp2ncd_CPPFLAGS}")
endif(USE_PARALLEL_COMPILER)

set_property(TARGET bp2ncd PROPERTY LINKER_LANGUAGE CXX)

#install(FILES bp2ncd.h DESTINATION ${PROJECT_BINARY_DIR}/utils/bp2ncd)
install(PROGRAMS ${PROJECT_BINARY_DIR}/utils/bp2ncd/bp2ncd DESTINATION ${bindir})

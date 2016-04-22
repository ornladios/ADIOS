include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/utils/adios_list_methods)
link_directories(${PROJECT_BINARY_DIR}/utils/adios_list_methods)

if(BUILD_WRITE)
  add_executable(adios_list_methods adios_list_methods.c)
  target_link_libraries(adios_list_methods adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES})
  set_target_properties(adios_list_methods PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}WRITE ${ADIOSLIB_CPPFLAGS} ${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${MPI_C_COMPILE_FLAGS}")
  set_property(TARGET adios_list_methods PROPERTY LINKER_LANGUAGE CXX)

  if(MPI_LINK_FLAGS)
     set_target_properties(adios_list_methods PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()

  install(PROGRAMS ${PROJECT_BINARY_DIR}/utils/adios_list_methods/adios_list_methods DESTINATION ${bindir})
endif(BUILD_WRITE)

add_executable(adios_list_methods_readonly adios_list_methods.c)
target_link_libraries(adios_list_methods_readonly adiosread ${ADIOSREADLIB_LDADD} ${MPI_C_LIBRARIES})
set_target_properties(adios_list_methods_readonly PROPERTIES COMPILE_FLAGS
"${ADIOSREADLIB_CPPFLAGS} ${ADIOSREADLIB_CFLAGS} ${ADIOSREADLIB_EXTRA_CPPFLAGS} ${MPI_C_COMPILE_FLAGS}")
set_property(TARGET adios_list_methods_readonly PROPERTY LINKER_LANGUAGE CXX)

if(MPI_LINK_FLAGS)
    set_target_properties(adios_list_methods_readonly PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()

install(PROGRAMS ${PROJECT_BINARY_DIR}/utils/adios_list_methods/adios_list_methods_readonly DESTINATION ${bindir})

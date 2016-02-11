include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/utils/list_methods)
link_directories(${PROJECT_BINARY_DIR}/utils/list_methods)

if(BUILD_WRITE)
  add_executable(list_methods list_methods.c)
  target_link_libraries(list_methods adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES})
  set_target_properties(list_methods PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}WRITE ${ADIOSLIB_CPPFLAGS} ${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${MPI_C_COMPILE_FLAGS}")
  set_property(TARGET list_methods PROPERTY LINKER_LANGUAGE CXX)


  if(MPI_LINK_FLAGS)
     set_target_properties(list_methods PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()

  install(PROGRAMS ${PROJECT_BINARY_DIR}/utils/list_methods/list_methods DESTINATION ${bindir})
endif(BUILD_WRITE)

add_executable(list_methods_readonly list_methods.c)
target_link_libraries(list_methods_readonly adiosread ${ADIOSREADLIB_LDADD} ${MPI_C_LIBRARIES})
set_target_properties(list_methods_readonly PROPERTIES COMPILE_FLAGS
"${ADIOSREADLIB_CPPFLAGS} ${ADIOSREADLIB_CFLAGS} ${ADIOSREADLIB_EXTRA_CPPFLAGS} ${MPI_C_COMPILE_FLAGS}")

if(MPI_LINK_FLAGS)
    set_target_properties(list_methods_readonly PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()

install(PROGRAMS ${PROJECT_BINARY_DIR}/utils/list_methods/list_methods_readonly DESTINATION ${bindir})

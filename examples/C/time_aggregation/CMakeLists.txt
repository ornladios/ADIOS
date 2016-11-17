include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/time_aggregation)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/time_aggregation)

add_executable(time_aggregation time_aggregation.c)
target_link_libraries(time_aggregation adios ${ADIOSLIB_LDADD})

add_executable(time_aggregation_no_xml time_aggregation_no_xml.c)
target_link_libraries(time_aggregation_no_xml adios ${ADIOSLIB_LDADD})

set (PROGS time_aggregation time_aggregation_no_xml)
foreach (PROG ${PROGS} )
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROG} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
  target_link_libraries(${PROG} adios ${MPI_C_LIBRARIES})
endforeach()

file(COPY time_aggregation.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/C/time_aggregation)

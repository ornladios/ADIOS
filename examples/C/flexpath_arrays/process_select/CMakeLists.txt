include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/flexpath_arrays/process_select)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/flexpath_arrays/process_select)

add_executable(select_arrays_write arrays_write.c)
add_executable(select_arrays_read arrays_read.c)

set (PROGS select_arrays_write select_arrays_read)
foreach (PROG ${PROGS} )
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROG} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
  target_link_libraries(${PROG} adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES})
endforeach()

file(COPY arrays.xml job.pbs 
     DESTINATION ${PROJECT_BINARY_DIR}/examples/C/flexpath_arrays/process_select)

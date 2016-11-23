include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/flexpath_arrays/process_select)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/flexpath_arrays/process_select)

add_executable(select_arrays_write arrays_write.c)
add_executable(select_arrays_read arrays_read.c)

set (PROGS_READ select_arrays_read)
set (PROGS_WRITE select_arrays_write)

foreach (PROG ${PROGS_WRITE})
  set_property(TARGET ${PROG} PROPERTY COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
  target_link_libraries(${PROG} adios ${ADIOSLIB_LDADD})
endforeach()

foreach (PROG ${PROGS_READ})
  set_property(TARGET ${PROG} PROPERTY COMPILE_FLAGS "${ADIOSREADLIB_CFLAGS}")
  target_link_libraries(${PROG} adios ${ADIOSREADLIB_LDADD})
endforeach()

file(COPY arrays.xml job.pbs
     DESTINATION ${PROJECT_BINARY_DIR}/examples/C/flexpath_arrays/process_select)

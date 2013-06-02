include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/arrays)
link_directories(${PROJECT_BINARY_DIR}/examples/C/arrays)

add_executable(arrays_write arrays_write.c)
target_link_libraries(arrays_write adios ${ADIOSLIB_LDADD})

add_executable(arrays_read arrays_read.c)
target_link_libraries(arrays_read adiosread ${ADIOSREADLIB_LDADD})

#install(FILES arrays.xml gwrite_arrays.ch gread_arrays.ch job.pbs  DESTINATION ${PROJECT_BINARY_DIR}/examples/C/arrays)


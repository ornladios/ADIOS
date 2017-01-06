include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/scalars)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/scalars)

add_executable(scalars_write scalars_write.c)
set_target_properties(scalars_write PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(scalars_write adios ${ADIOSLIB_LDADD})

add_executable(scalars_read scalars_read.c)
set_target_properties(scalars_read PROPERTIES COMPILE_FLAGS "${ADIOSREADLIB_CFLAGS}")
target_link_libraries(scalars_read adiosread ${ADIOSREADLIB_LDADD})

add_executable(scalars_write_step scalars_write.c)
set_target_properties(scalars_write_step PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(scalars_write_step adios ${ADIOSLIB_LDADD})

#install(FILES scalars.xml gwrite_scalars.ch gread_scalars.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/scalars)

include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/scalars)
link_directories(${PROJECT_BINARY_DIR}/examples/C/scalars)

add_executable(scalars_write scalars_write.c)
target_link_libraries(scalars_write ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(scalars_read scalars_read.c)
target_link_libraries(scalars_read ${ADIOS_LIB} ${ADIOSLIB_LDADD})

install(FILES scalars.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/C/scalars)

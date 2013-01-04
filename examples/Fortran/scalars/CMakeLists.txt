include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/scalars)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/scalars)

add_executable(scalars_write_f scalars_write.F90)
target_link_libraries(scalars_write_f ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(scalars_read_f scalars_read.F90)
target_link_libraries(scalars_read_f ${ADIOS_LIB} ${ADIOSLIB_LDADD})

install(FILES scalars.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/scalars)

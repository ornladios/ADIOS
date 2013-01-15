include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/arrays)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/arrays)

add_executable(arrays_write_f arrays_write.F90)
target_link_libraries(arrays_write_f adiosf ${ADIOSLIB_LDADD})

add_executable(arrays_read_f arrays_read.F90)
target_link_libraries(arrays_read_f adiosf ${ADIOSLIB_LDADD})

install(FILES arrays.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/arrays)

include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/arrays)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/arrays)

add_executable(arrays_write arrays_write.F90)
target_link_libraries(arrays_write ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(arrays_read arrays_read.F90)
target_link_libraries(arrays_read ${ADIOS_LIB} ${ADIOSLIB_LDADD})

install(FILES arrays.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/arrays)

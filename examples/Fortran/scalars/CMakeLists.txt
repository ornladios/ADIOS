include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/scalars)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/scalars)

add_executable(scalars_write_f scalars_write.F90)
target_link_libraries(scalars_write_f adiosf ${ADIOSLIB_LDADD})

add_executable(scalars_read_f scalars_read.F90)
target_link_libraries(scalars_read_f adiosreadf ${ADIOSREADLIB_LDADD})

#install(FILES scalars.xml gwrite_scalars.fh gread_scalars.fh job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/scalars)

include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/global-array-time)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/global-array-time)

add_executable(adios_globaltime_f adios_globaltime.F90)
target_link_libraries(adios_globaltime_f adiosf ${ADIOSLIB_LDADD})

#install(FILES adios_globaltime.xml gread_restart.fh gwrite_restart.fh job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/global-array-time)

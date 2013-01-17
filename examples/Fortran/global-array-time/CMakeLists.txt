include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/global-array-time)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/global-array-time)

add_executable(adios_globaltime_f adios_globaltime.F90)
target_link_libraries(adios_globaltime_f adiosf ${ADIOSLIB_LDADD})

install(FILES global-array-time.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/global-array-time)

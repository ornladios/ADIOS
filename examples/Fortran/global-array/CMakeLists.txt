include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/global-array)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/global-array)

add_executable(adios_global adios_global.F90)
target_link_libraries(adios_global ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(adios_global_no_xml adios_global_no_xml.F90)
target_link_libraries(adios_global_no_xml ${ADIOS_LIB} ${ADIOSLIB_LDADD})

install(FILES adios_global.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/global-array)

include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src  ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/global-array)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/global-array)

add_executable(adios_global_f adios_global.F90)
target_link_libraries(adios_global_f adiosf ${ADIOSLIB_LDADD})

add_executable(adios_global_no_xml_f adios_global_no_xml.F90)
target_link_libraries(adios_global_no_xml_f adiosf ${ADIOSLIB_LDADD})

#install(FILES adios_global.xml gread_temperature.fh gwrite_temperature.fh job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/global-array)

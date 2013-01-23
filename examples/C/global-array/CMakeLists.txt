include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/global-array)
link_directories(${PROJECT_BINARY_DIR}/examples/C/global-array)

add_executable(adios_global adios_global.c)
target_link_libraries(adios_global adios ${ADIOSLIB_LDADD})

add_executable(adios_global_2files adios_global_2files.c)
target_link_libraries(adios_global_2files adios ${ADIOSLIB_LDADD})

add_executable(adios_read_global adios_read_global.c)
target_link_libraries(adios_read_global adiosread ${ADIOSLIB_LDADD})
set_target_properties(adios_read_global PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

add_executable(adios_read_writeblock adios_read_writeblock.c)
target_link_libraries(adios_read_writeblock adiosread ${ADIOSLIB_LDADD})

add_executable(adios_read_chunk adios_read_chunk.c)
target_link_libraries(adios_read_chunk adiosread ${ADIOSLIB_LDADD})

add_executable(adios_read_gpp adios_read_gpp.c)
target_link_libraries(adios_read_gpp adiosread ${ADIOSLIB_LDADD})

add_executable(adios_global_no_xml adios_global_no_xml.c)
target_link_libraries(adios_global_no_xml adios ${ADIOSLIB_LDADD})

add_executable(adios_read_global_no_xml adios_read_global_no_xml.c)
target_link_libraries(adios_read_global_no_xml adiosread ${ADIOSLIB_LDADD})
set_target_properties(adios_read_global_no_xml PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

install(FILES adios_global.xml gread_temperature.ch gwrite_temperature.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/global-array)

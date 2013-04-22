include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/global-array-time)
link_directories(${PROJECT_BINARY_DIR}/examples/C/global-array-time)

add_executable(adios_globaltime adios_globaltime.c)
target_link_libraries(adios_globaltime adios ${ADIOSLIB_LDADD})

add_executable(adios_globaltime_no_xml adios_globaltime_no_xml.c)
target_link_libraries(adios_globaltime_no_xml adios ${ADIOSLIB_LDADD})

add_executable(adios_read_globaltime adios_read_globaltime.c)
target_link_libraries(adios_read_globaltime adiosread ${ADIOSLIB_LDADD})
set_target_properties(adios_read_globaltime PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

add_executable(adios_read_globaltime_streaming adios_read_globaltime_streaming.c)
target_link_libraries(adios_read_globaltime_streaming adiosread ${ADIOSLIB_LDADD})

#install(FILES adios_globaltime.xml gread_restart.ch gwrite_restart.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/global-array-time)

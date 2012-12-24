include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/global-array)
link_directories(${PROJECT_BINARY_DIR}/examples/C/global-array)

add_executable(adios_global adios_global.c)
target_link_libraries(adios_global ${ADIOS_LIB} ${ADIOSLIB_LDADD} ${MPI_LIBRARIES})

add_executable(adios_global_2files adios_global_2files.c)
target_link_libraries(adios_global_2files ${ADIOS_LIB} ${ADIOSLIB_LDADD} ${MPI_LIBRARIES})

add_executable(adios_read_global adios_read_global.c)
target_link_libraries(adios_read_global ${ADIOSREAD_LIB} ${ADIOSLIB_LDADD} ${MPI_LIBRARIES})
set_target_properties(adios_read_global PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

add_executable(adios_read_writeblock adios_read_writeblock.c)
target_link_libraries(adios_read_writeblock ${ADIOSREAD_LIB} ${ADIOSLIB_LDADD}
${MPI_LIBRARIES})

add_executable(adios_read_chunk adios_read_chunk.c)
target_link_libraries(adios_read_chunk ${ADIOSREAD_LIB} ${ADIOSLIB_LDADD}
${MPI_LIBRARIES})

add_executable(adios_read_gpp adios_read_gpp.c)
target_link_libraries(adios_read_gpp ${ADIOSREAD_LIB} ${ADIOSLIB_LDADD}
${MPI_LIBRARIES})

add_executable(adios_global_no_xml adios_global_no_xml.c)
target_link_libraries(adios_global_no_xml ${ADIOS_LIB} ${ADIOSLIB_LDADD}
${MPI_LIBRARIES})

add_executable(adios_read_global_no_xml adios_read_global_no_xml.c)
target_link_libraries(adios_read_global_no_xml ${ADIOSREAD_LIB}
${ADIOSLIB_LDADD} ${MPI_LIBRARIES})
set_target_properties(adios_read_global_no_xml PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

install(FILES adios_global.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/C/global-array)

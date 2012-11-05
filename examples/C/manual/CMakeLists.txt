include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/manual)
link_directories(${PROJECT_BINARY_DIR}/examples/C/manual)

add_executable(1_nonadios_example 1_nonadios_example.c)
target_link_libraries(1_nonadios_example ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(2_adios_write 2_adios_write.c)
target_link_libraries(2_adios_write ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(3_adios_read 3_adios_read.c)
target_link_libraries(3_adios_read ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(4_adios_nfiles 4_adios_nfiles.c)
target_link_libraries(4_adios_nfiles ${ADIOS_LIB} ${ADIOSLIB_LDADD})

install(FILES config.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/C/manual)


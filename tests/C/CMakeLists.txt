include_directories(${PROJECT_SOURCE_DIR}/tests/C)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${MXML_INCLUDE_DIR}/../lib)

add_executable(adios_test_c adios_test_c.c)
target_link_libraries(adios_test_c adios ${ADIOSLIB_LDADD}) 
set(adios_test_c_CPPFLAGS "${ADIOSLIB_CPPFLAGS}")
set(adios_test_c_CFLAGS "${ADIOSLIB_CFLAGS}")
set_target_properties(adios_test_c PROPERTIES COMPILE_FLAGS "${adios_test_c_CFLAGS} ${adios_test_c_CPPFLAGS}")

install(FILES config_c.xml DESTINATION ${PROJECT_BINARY_DIR}/tests/C)

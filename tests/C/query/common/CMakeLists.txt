include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)
include_directories(${MXML_INCLUDE_DIR})
link_directories(${PROJECT_BINARY_DIR}/src)

add_library(parse_test_query_xml adios_query_xml_parse.c adios_query_xml_parse.h)
set_target_properties(parse_test_query_xml PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(parse_test_query_xml adios)

add_executable(compute_expected_query_results compute_expected_query_results.c)
set_target_properties(compute_expected_query_results PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(compute_expected_query_results parse_test_query_xml adios ${ADIOSLIB_LDADD}) 

add_executable(adios_query_test adios_query_test.c)
target_link_libraries(adios_query_test parse_test_query_xml adios ${ADIOSLIB_LDADD})
set_target_properties(adios_query_test PROPERTIES COMPILE_FLAGS "${ALACRITY_CPPFLAGS} ${ADIOSLIB_CFLAGS}")

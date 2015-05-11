include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/query)

add_executable(write_vars write_vars.c)
target_link_libraries(write_vars  adios_nompi ${ADIOSLIB_SEQ_LDADD})
set_target_properties(write_vars  PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}_NOMPI")

add_executable(query_vars query_vars.c)
target_link_libraries(query_vars  adiosread_nompi ${ADIOSREADLIB_SEQ_LDADD})
set_target_properties(query_vars  PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}_NOMPI")

add_executable(write_table write_table.c)
target_link_libraries(write_table adios_nompi ${ADIOSLIB_SEQ_LDADD})
set_target_properties(write_table PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}_NOMPI")

add_executable(query_table query_table.c)
target_link_libraries(query_table adiosread_nompi ${ADIOSREADLIB_SEQ_LDADD})
set_target_properties(query_table PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}_NOMPI")


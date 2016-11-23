include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test/C/flexpath_tests/include)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)


add_executable(arrays_writer writer.c)
set_property(TARGET arrays_writer PROPERTY COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(arrays_writer adios ${ADIOSLIB_LDADD})

add_executable(arrays_reader reader.c)
set_property(TARGET arrays_reader PROPERTY COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(arrays_reader adios ${ADIOSLIB_LDADD})

file(COPY arrays.xml DESTINATION ${PROJECT_BINARY_DIR}/tests/C/flexpath_tests/global_range_select)

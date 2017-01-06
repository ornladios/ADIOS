include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/tests/C/flexpath_tests/include)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)


add_executable(writer_maya_noxml writer.c ../common/utils.c)
set_property(TARGET writer_maya_noxml PROPERTY COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(writer_maya_noxml adios ${ADIOSLIB_LDADD})

add_executable(reader_maya_noxml reader.c ../common/utils.c)
set_property(TARGET reader_maya_noxml PROPERTY COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(reader_maya_noxml adios ${ADIOSLIB_LDADD})

file(COPY readme.txt DESTINATION ${PROJECT_BINARY_DIR}/tests/C/flexpath_tests/maya_append)

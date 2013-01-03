include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/read_all)
link_directories(${PROJECT_BINARY_DIR}/examples/C/read_all)

add_executable(read_all read_all.c)
target_link_libraries(read_all ${ADIOSREAD_NOMPI_LIB} ${ADIOSLIB_LDADD})
set_target_properties(read_all PROPERTIES COMPILE_FLAGS "-D_NOMPI -DADIOS_USE_READ_API_1")

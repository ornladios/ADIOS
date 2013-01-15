include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/attributes)
link_directories(${PROJECT_BINARY_DIR}/examples/C/attributes)

add_executable(attributes_write attributes_write.c)
target_link_libraries(attributes_write adios ${ADIOSLIB_LDADD})

add_executable(attributes_read attributes_read.c)
target_link_libraries(attributes_read adiosread ${ADIOSLIB_LDADD})
set_target_properties(attributes_read PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

install(FILES attributes.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/C/attributes)


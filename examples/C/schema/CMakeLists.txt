include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/schema)
link_directories(${PROJECT_BINARY_DIR}/examples/C/schema)

add_executable(uniform uniform.c)
target_link_libraries(uniform ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(rectilinear rectilinear.c)
target_link_libraries(rectilinear ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(structured structured.c)
target_link_libraries(structured ${ADIOS_LIB} ${ADIOSLIB_LDADD})

add_executable(unstructured unstructured.c)
target_link_libraries(unstructured ${ADIOS_LIB} ${ADIOSLIB_LDADD})

install(FILES structured.xml rectilinear.xml unstructured.xml uniform.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/C/schema)

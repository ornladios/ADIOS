include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(bpdump bpdump.c)
target_link_libraries(bpdump ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD})

install(FILES bpdump DESTINATION
${PROJECT_BINARY_DIR}/bin)

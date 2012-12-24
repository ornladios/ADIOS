include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(bp2ascii bp2ascii.c)
target_link_libraries(bp2ascii ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD} "pthread")

install(FILES bp2ascii DESTINATION
${PROJECT_BINARY_DIR}/bin)

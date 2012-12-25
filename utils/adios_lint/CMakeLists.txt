include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(adios_lint adios_lint.c)
target_link_libraries(adios_lint ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD} "pthread")

install(FILES adios_lint DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

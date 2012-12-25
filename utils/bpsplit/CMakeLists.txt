include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(bpsplit bpsplit.c)
target_link_libraries(bpsplit ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD} "pthread")

add_executable(bpappend bpappend.c)
target_link_libraries(bpappend ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD} "pthread")

add_executable(bpgettime bpgettime.c)
target_link_libraries(bpgettime ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD} "pthread")

install(FILES bpsplit bpappend bpgettime DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/utils/bpsplit)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
link_directories(${PROJECT_BINARY_DIR}/utils/bpsplit)

add_executable(bpsplit bpsplit.c)
target_link_libraries(bpsplit ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD})
set_target_properties(bpsplit PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

add_executable(bpappend bpappend.c)
target_link_libraries(bpappend  ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD})
set_target_properties(bpappend PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

add_executable(bpgettime bpgettime.c)
target_link_libraries(bpgettime ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD})
set_target_properties(bpgettime PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

if(USE_PARALLEL_COMPILER)
endif(USE_PARALLEL_COMPILER)

install(FILES bpsplit bpappend bpgettime DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

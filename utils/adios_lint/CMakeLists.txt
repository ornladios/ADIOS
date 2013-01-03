include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/utils/adios_lint)
link_directories(${PROJECT_BINARY_DIR}/utils/adios_lint)

add_executable(adios_lint adios_lint.c)
target_link_libraries(adios_lint ${ADIOS_INTERNAL_NOMPI_LIB} ${ADIOSLIB_INT_LDADD})
set_target_properties(adios_lint PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

#install(FILES config.xml DESTINATION ${PROJECT_BINARY_DIR}/utils/adios_lint)

if(USE_PARALLEL_COMPILER)
endif(USE_PARALLEL_COMPILER)

install(FILES adios_lint DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

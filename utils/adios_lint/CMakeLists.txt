include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/utils/adios_lint)
include_directories(${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/utils/adios_lint)

add_executable(adios_lint adios_lint.c)
target_link_libraries(adios_lint adios_internal_nompi ${ADIOSLIB_INT_LDADD})
set_target_properties(adios_lint PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")
set_property(TARGET adios_lint PROPERTY LINKER_LANGUAGE CXX)

#install(FILES config.xml DESTINATION ${PROJECT_BINARY_DIR}/utils/adios_lint)
install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/adios_lint/adios_lint DESTINATION ${bindir})

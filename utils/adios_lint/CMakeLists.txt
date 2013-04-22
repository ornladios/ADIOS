include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/utils/adios_lint)
link_directories(${PROJECT_BINARY_DIR}/utils/adios_lint)

add_executable(adios_lint adios_lint.c)
target_link_libraries(adios_lint adios_internal_nompi ${ADIOSLIB_INT_LDADD})
set_target_properties(adios_lint PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

#install(FILES config.xml DESTINATION ${PROJECT_BINARY_DIR}/utils/adios_lint)
install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/adios_lint/adios_lint DESTINATION ${bindir})

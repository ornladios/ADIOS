include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/core/transforms)
include_directories(${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(bpmeta bpmeta.c)
target_link_libraries(bpmeta adios_internal_nompi ${ADIOSLIB_INT_LDADD})
set_target_properties(bpmeta PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")
set_property(TARGET bpmeta PROPERTY LINKER_LANGUAGE CXX)

install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/bpmeta/bpmeta DESTINATION ${bindir})

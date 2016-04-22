include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/core/transforms)
include_directories(${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(bpdump bpdump.c)
target_link_libraries(bpdump adios_internal_nompi ${ADIOSLIB_INT_LDADD})
set_target_properties(bpdump PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")
set_property(TARGET bpdump PROPERTY LINKER_LANGUAGE CXX)

install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/bpdump/bpdump DESTINATION ${bindir})

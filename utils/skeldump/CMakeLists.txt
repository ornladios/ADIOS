include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_BINARY_DIR}  ${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src/core)

add_executable(skeldump skeldump.c)
#target_link_libraries(skeldump adios_internal_nompi ${ADIOSLIB_INT_LDADD})
target_link_libraries(skeldump adiosread_nompi ${ADIOSREADLIB_SEQ_LDADD})
set_target_properties(skeldump PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/skeldump/skeldump DESTINATION ${bindir})

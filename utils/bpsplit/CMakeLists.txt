include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/utils/bpsplit)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
link_directories(${PROJECT_BINARY_DIR}/utils/bpsplit)

add_executable(bpsplit bpsplit.c)
target_link_libraries(bpsplit adios_internal_nompi ${ADIOSLIB_INT_LDADD})
set_target_properties(bpsplit PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

add_executable(bpappend bpappend.c)
target_link_libraries(bpappend  adios_internal_nompi ${ADIOSLIB_INT_LDADD})
set_target_properties(bpappend PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

add_executable(bpgettime bpgettime.c)
target_link_libraries(bpgettime adios_internal_nompi ${ADIOSLIB_INT_LDADD})
set_target_properties(bpgettime PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")


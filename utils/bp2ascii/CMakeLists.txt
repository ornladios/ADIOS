include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/utils/bp2ascii)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
link_directories(${PROJECT_BINARY_DIR}/utils/bp2ascii)

add_executable(bp2ascii bp2ascii.c)
target_link_libraries(bp2ascii adios_internal_nompi ${ADIOSLIB_INT_LDADD})
set_target_properties(bp2ascii PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_INT_CFLAGS}")

if(USE_PARALLEL_COMPILER)
endif(USE_PARALLEL_COMPILER)                             

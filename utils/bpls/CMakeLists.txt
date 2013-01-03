include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/utils/bpls)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
link_directories(${PROJECT_BINARY_DIR}/utils/bpls)

add_executable(bpls bpls.c)
target_link_libraries(bpls ${ADIOSREAD_NOMPI_LIB} ${ADIOSREADLIB_SEQ_LDADD})
set_target_properties(bpls PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_SEQ_CPPFLAGS} ${ADIOSREADLIB_SEQ_CFLAGS}")

if(USE_PARALLEL_COMPILER)
endif(USE_PARALLEL_COMPILER)

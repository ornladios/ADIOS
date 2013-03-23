include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/utils/bp2bp)
link_directories(${PROJECT_BINARY_DIR}/utils/bp2bp)

add_executable(bp2bp bp2bp.c)
target_link_libraries(bp2bp adios ${ADIOSLIB_LDADD})
set_target_properties(bp2bp PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1 ${ADIOSLIB_CPPFLAGS} ${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS}")

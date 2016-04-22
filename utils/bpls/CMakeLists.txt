include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/utils/bpls)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/utils/bpls)

add_executable(bpls bpls.c)
target_link_libraries(bpls adiosread_nompi ${ADIOSREADLIB_SEQ_LDADD})
set_target_properties(bpls PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_SEQ_CPPFLAGS} ${ADIOSREADLIB_SEQ_CFLAGS}")
set_property(TARGET bpls PROPERTY LINKER_LANGUAGE CXX)

#install(FILES bpls.h DESTINATION ${PROJECT_BINARY_DIR}/utils/bpls)
install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/bpls/bpls DESTINATION ${bindir})

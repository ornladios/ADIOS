include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/utils/bpdiff)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/utils/bpdiff)

add_executable(bpdiff bpdiff.c decompose_block.c utils.c)
target_link_libraries(bpdiff adiosread_nompi ${ADIOSREADLIB_LDADD})
set_target_properties(bpdiff PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS} ${ADIOSREADLIB_CFLAGS}")
set_property(TARGET bpdiff PROPERTY LINKER_LANGUAGE CXX)

install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/bpdiff/bpdiff DESTINATION ${bindir})

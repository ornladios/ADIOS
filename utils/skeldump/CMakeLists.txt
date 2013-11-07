include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/utils/skeldump)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
link_directories(${PROJECT_BINARY_DIR}/utils/skeldump)

add_executable(skeldump skeldump.c)
target_link_libraries(skeldump adiosread_nompi ${ADIOSREADLIB_SEQ_LDADD})
set_target_properties(skeldump PROPERTIES COMPILE_FLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_SEQ_CPPFLAGS} ${ADIOSREADLIB_SEQ_CFLAGS}")

#install(FILES skeldump.h DESTINATION ${PROJECT_BINARY_DIR}/utils/skeldump)
install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/skeldump/skeldump DESTINATION ${bindir})

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/core/transforms)
include_directories(${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(adios_index_fastbit adios_index_fastbit.c)
target_link_libraries(adios_index_fastbit adios ${ADIOSLIB_LDADD})
set_target_properties(adios_index_fastbit PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CPPFLAGS} ${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS}")

install(PROGRAMS ${CMAKE_BINARY_DIR}/utils/fastbit/adios_index_fastbit DESTINATION ${bindir})

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(bp2bp bp2bp.c)
target_link_libraries(bp2bp ${ADIOS_LIB} ${ADIOSLIB_LDADD} ${MPI_LIBRARIES})
set_target_properties(bp2bp PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1 ${ADIOSLIB_CPPFLAGS} ${ADIOSLIB_CFLAGS} ${CPPFLAGS}")

install(FILES bp2bp DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

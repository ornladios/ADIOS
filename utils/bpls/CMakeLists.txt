include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(bpls bpls.c)
target_link_libraries(bpls ${ADIOSREAD_NOMPI_LIB} ${ADIOSREADLIB_SEQ_LDADD})
set_target_properties(bpls PROPERTIES COMPILE_FLAGS ${ADIOSREADLIB_SEQ_CPPFLAGS})

install(FILES bpls DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


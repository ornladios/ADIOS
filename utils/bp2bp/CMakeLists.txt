include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/utils/bp2bp)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/utils/bp2bp)

add_executable(bp2bp bp2bp.c)
target_link_libraries(bp2bp adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES})
set_target_properties(bp2bp PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}ADIOS_USE_READ_API_1 ${ADIOSLIB_CPPFLAGS} ${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${MPI_C_COMPILE_FLAGS}")
set_property(TARGET bp2bp PROPERTY LINKER_LANGUAGE CXX)

if(MPI_LINK_FLAGS)
   set_target_properties(bp2bp PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()

install(PROGRAMS ${PROJECT_BINARY_DIR}/utils/bp2bp/bp2bp DESTINATION ${bindir})

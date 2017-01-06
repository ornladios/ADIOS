include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)


add_executable(posix_fgr posix_fgr.c)
set_target_properties(posix_fgr PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
target_link_libraries(posix_fgr adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES}) 


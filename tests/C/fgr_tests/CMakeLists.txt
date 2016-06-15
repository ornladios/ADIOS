include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)
link_directories(${MXML_LIB_DIR})


add_executable(posix_fgr posix_fgr.c)
target_link_libraries(posix_fgr adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES}) 


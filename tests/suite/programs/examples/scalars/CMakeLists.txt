include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/scalars)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/test/suite/programs/examples/scalars)

add_executable(scalars_write_C scalars_write_C.c)
if(MPI_COMPILE_FLAGS)
  set_target_properties(scalars_write_C PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(scalars_write_C PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()
target_link_libraries(scalars_write_C adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES})

add_executable(scalars_read_C scalars_read_C.c)
if(MPI_COMPILE_FLAGS)
  set_target_properties(scalars_read_C PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(scalars_read_C PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()
target_link_libraries(scalars_read_C adiosread ${ADIOSREADLIB_LDADD} ${MPI_C_LIBRARIES})

if(BUILD_FORTRAN)
    add_executable(scalars_write_F scalars_write_F.F90)
    target_link_libraries(scalars_write_F adiosf ${ADIOSLIB_LDADD})

    add_executable(scalars_read_F scalars_read_F.F90)
    target_link_libraries(scalars_read_F adiosreadf ${ADIOSREADLIB_LDADD})
endif

#install(FILES scalars.xml gwrite_scalars.ch gread_scalars.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/scalars)

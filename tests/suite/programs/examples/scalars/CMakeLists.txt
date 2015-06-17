include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test/suite/programs/examples/scalars)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/test/suite/programs/examples/scalars)

add_executable(scalars_write_C scalars_write_C.c gwrite_scalars.ch)
if(MPI_COMPILE_FLAGS)
  set_target_properties(scalars_write_C PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(scalars_write_C PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()
target_link_libraries(scalars_write_C adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES})
add_custom_command(
	OUTPUT gwrite_scalars.ch
	COMMAND ${PROJECT_SOURCE_DIR}/utils/gpp/gpp.py ${PROJECT_SOURCE_DIR}/tests/suite/programs/examples/scalars/scalars_C.xml
	DEPENDS scalars_C.xml
	)

add_executable(scalars_read_C scalars_read_C.c)
if(MPI_COMPILE_FLAGS)
  set_target_properties(scalars_read_C PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(scalars_read_C PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()
target_link_libraries(scalars_read_C adiosread ${ADIOSREADLIB_LDADD} ${MPI_C_LIBRARIES})

if(BUILD_FORTRAN)
    add_executable(scalars_write_F scalars_write_F.F90 gwrite_scalars.fh)
    target_link_libraries(scalars_write_F adiosf ${ADIOSLIB_LDADD})
    add_custom_command(
	OUTPUT gwrite_scalars.fh
	COMMAND ${PROJECT_SOURCE_DIR}/utils/gpp/gpp.py ${PROJECT_SOURCE_DIR}/tests/suite/programs/examples/scalars/scalars_F.xml
	DEPENDS scalars_F.xml
	)


    add_executable(scalars_read_F scalars_read_F.F90)
    target_link_libraries(scalars_read_F adiosreadf ${ADIOSREADLIB_LDADD})
endif

#install(FILES scalars.xml gwrite_scalars.ch gread_scalars.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/scalars)

include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test/suite/programs/examples/global_array_time)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/test/suite/programs/examples/global_array_time)

add_executable(global_array_time_write_C global_array_time_write_C.c)
target_link_libraries(global_array_time_write_C adios ${ADIOSLIB_LDADD})
add_custom_command(
        OUTPUT gwrite_restart.ch
        COMMAND ${PROJECT_SOURCE_DIR}/utils/gpp/gpp.py
        ${PROJECT_SOURCE_DIR}/tests/suite/programs/examples/global_array_time/global_array_time_C.xml
        DEPENDS global_array_C.xml
        )

add_executable(global_array_time_read_as_file_C global_array_time_read_as_file_C.c)
target_link_libraries(global_array_time_read_as_file_C adiosread ${ADIOSREADLIB_LDADD})
#set_target_properties(global_array_time_read_as_file_C PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

add_executable(global_array_time_read_as_stream_C global_array_time_read_as_stream_C.c)
target_link_libraries(global_array_time_read_as_stream_C adiosread ${ADIOSREADLIB_LDADD})

set (PROGS global_array_time_read_as_stream_C global_array_time_read_as_file_C global_array_time_write_C)
foreach (PROG ${PROGS} )
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROG} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
  target_link_libraries(${PROG} adios ${MPI_C_LIBRARIES})
endforeach()


if(BUILD_FORTRAN)
    add_executable(global_array_time_write_F global_array_time_write_F.F90)
    target_link_libraries(global_array_time_write_F adiosf ${ADIOSLIB_LDADD})
    add_custom_command(
        OUTPUT gwrite_restart.fh
        COMMAND ${PROJECT_SOURCE_DIR}/utils/gpp/gpp.py
        ${PROJECT_SOURCE_DIR}/tests/suite/programs/examples/global_array_time/global_array_time_F.xml
        DEPENDS global_array_F.xml
        )
endif

#install(FILES global_array_time_write_C.xml gread_restart.ch gwrite_restart.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/global-array-time)

include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test/suite/programs/examples/global_array)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/test/suite/programs/examples/global_array)

add_executable(global_array_write_C global_array_write_C.c)
target_link_libraries(global_array_write_C adios ${ADIOSLIB_LDADD})
add_custom_command(
	OUTPUT gwrite_temperature.ch
	COMMAND ${PROJECT_SOURCE_DIR}/utils/gpp/gpp.py ${PROJECT_SOURCE_DIR}/tests/suite/programs/examples/global_array/global_array_C.xml
	DEPENDS global_array_C.xml
	)

add_executable(global_array_write_noxml_C global_array_write_noxml_C.c)
target_link_libraries(global_array_write_noxml_C adios ${ADIOSLIB_LDADD})

add_executable(global_array_write_byid_noxml_C global_array_write_byid_noxml_C.c)
target_link_libraries(global_array_write_byid_noxml_C adios ${ADIOSLIB_LDADD})

add_executable(global_array_aggregate_by_color_C global_array_aggregate_by_color_C.c)
target_link_libraries(global_array_aggregate_by_color_C adios ${ADIOSLIB_LDADD})

add_executable(global_array_read_C global_array_read_C.c)
target_link_libraries(global_array_read_C adiosread ${ADIOSREADLIB_LDADD})

add_executable(global_array_read_noxml_C global_array_read_noxml_C.c)
target_link_libraries(global_array_read_noxml_C adiosread ${ADIOSREADLIB_LDADD})

add_executable(global_array_read_byid_noxml_C global_array_read_byid_noxml_C.c)
target_link_libraries(global_array_read_byid_noxml_C adiosread ${ADIOSREADLIB_LDADD})

#set_target_properties(global_array_read_noxml_C PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

set (PROGS global_array_write_C global_array_read_C global_array_write_noxml_C global_array_read_noxml_C)

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
    add_executable(global_array_write_F global_array_write_F.F90)
    target_link_libraries(global_array_write_F adiosf ${ADIOSLIB_LDADD})
    add_custom_command(
	OUTPUT gwrite_temperature.fh
	COMMAND ${PROJECT_SOURCE_DIR}/utils/gpp/gpp.py ${PROJECT_SOURCE_DIR}/tests/suite/programs/examples/global_array/global_array_F.xml
	DEPENDS global_array_F.xml
	)
    add_executable(global_array_write_noxml_F global_array_write_noxml_F.F90)
    target_link_libraries(global_array_write_noxml_F adiosf ${ADIOSLIB_LDADD})
    add_executable(global_array_write_byid_noxml_F global_array_write_byid_noxml_F.F90)
    target_link_libraries(global_array_write_byid_noxml_F adiosf ${ADIOSLIB_LDADD})
endif


#install(FILES adios_global.xml gread_temperature.ch gwrite_temperature.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/global-array)


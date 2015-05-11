include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/manual)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/manual)

set (PROGS 1_nonadios_example 2_adios_write 3_adios_read 4_adios_nfiles)

foreach (PROG ${PROGS} )
  add_executable(${PROG} ${PROG}.c)
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROG} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
  target_link_libraries(${PROG} adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES})
endforeach()

#install(FILES config.xml gread_temperature.ch gwrite_temperature.ch DESTINATION ${PROJECT_BINARY_DIR}/examples/C/manual)


include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/manual)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/manual)

set (PROGS 1_nonadios_example 2_adios_write 3_adios_read 4_adios_nfiles)

foreach (PROG ${PROGS} )
  add_executable(${PROG} ${PROG}.c)
  set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
  target_link_libraries(${PROG} adios ${ADIOSLIB_LDADD})
endforeach()

#install(FILES config.xml gread_temperature.ch gwrite_temperature.ch DESTINATION ${PROJECT_BINARY_DIR}/examples/C/manual)


include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/stat)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/stat)

add_executable(stat_write stat_write.c)
target_link_libraries(stat_write adios ${ADIOSLIB_LDADD})

add_executable(stat_read stat_read.c)
target_link_libraries(stat_read adios ${ADIOSREADLIB_LDADD})
#set_target_properties(stat_read PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

set (PROGS_READ stat_read)
set (PROGS_WRITE stat_write)

foreach (PROG ${PROGS_WRITE})
  set_property(TARGET ${PROG} PROPERTY COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
    target_link_libraries(${PROG} adios)
endforeach()

foreach (PROG ${PROGS_READ})
  set_property(TARGET ${PROG} PROPERTY COMPILE_FLAGS "${ADIOSREADLIB_CFLAGS}")
  target_link_libraries(${PROG} adios)
endforeach()

#install(FILES stat.xml gwrite_stat.ch gread_stat.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/stat)

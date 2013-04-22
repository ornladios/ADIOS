include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/stat)
link_directories(${PROJECT_BINARY_DIR}/examples/C/stat)

add_executable(stat_write stat_write.c)
target_link_libraries(stat_write adios ${ADIOSLIB_LDADD})

add_executable(stat_read stat_read.c)
target_link_libraries(stat_read adios ${ADIOSLIB_LDADD})
set_target_properties(stat_read PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

#install(FILES stat.xml gwrite_stat.ch gread_stat.ch job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/C/stat)

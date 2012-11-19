include_directories(${PROJECT_SOURCE_DIR}/tests/suite)
include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/tests/suite)


if(BUILD_WRITE)
  add_executable(local_array_time local_array_time.c)
  target_link_libraries(local_array_time ${ADIOS_LIB} ${ADIOSLIB_LDADD})

  add_executable(adios_amr_write adios_amr_write.c)
  target_link_libraries(adios_amr_write ${ADIOS_LIB} ${ADIOSLIB_LDADD})

  add_executable(adios_amr_write_2vars adios_amr_write_2vars.c)
  target_link_libraries(adios_amr_write_2vars ${ADIOS_LIB} ${ADIOSLIB_LDADD})
  
  add_executable(adios_staged_read adios_staged_read.c)
  target_link_libraries(adios_staged_read ${ADIOSREAD_LIB} ${ADIOSREADLIB_LDADD})

  add_executable(adios_staged_read_2vars adios_staged_read_2vars.c)
  target_link_libraries(adios_staged_read_2vars ${ADIOSREAD_LIB} ${ADIOSREADLIB_LDADD})

  add_executable(write_read write_read.c)
  target_link_libraries(write_read ${ADIOS_LIB} ${ADIOSLIB_LDADD})

  add_executable(write_alternate write_alternate.c)
  target_link_libraries(write_alternate ${ADIOS_LIB} ${ADIOSLIB_LDADD})

  add_executable(many_vars many_vars.c)
  target_link_libraries(many_vars ${ADIOS_LIB} ${ADIOSLIB_LDADD})

  add_executable(selections selections.c)
  target_link_libraries(selections ${ADIOS_LIB} ${ADIOSLIB_LDADD})
endif(BUILD_WRITE)

if(BUILD_FORTRAN)
  if(BUILD_WRITE)
    add_executable(posix_method posix_method.F90)
    target_link_libraries(posix_method ${ADIOS_LIB} ${ADIOSLIB_LDADD})
    #not complete
  endif(BUILD_WRITE)
endif(BUILD_FORTRAN)

install(FILES adios_amr_write_2vars.xml adios_amr_write.xml
local_array_time.xml posix_method.xml write_alternate.xml write_read.xml DESTINATION ${PROJECT_BINARY_DIR}/tests/suite/programs)

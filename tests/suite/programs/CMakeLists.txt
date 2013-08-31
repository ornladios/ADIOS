include_directories(${PROJECT_SOURCE_DIR}/tests/suite)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/tests/suite/programs)
include_directories(${PROJECT_BINARY_DIR}/tests/suite/programs)
link_directories(${PROJECT_BINARY_DIR}/tests/suite/programs)


if(BUILD_WRITE)
  add_executable(local_array_time local_array_time.c)
  target_link_libraries(local_array_time adios ${ADIOSLIB_LDADD})
  
  add_executable(adios_amr_write adios_amr_write.c)
  target_link_libraries(adios_amr_write adios ${ADIOSLIB_LDADD})

  add_executable(adios_amr_write_2vars adios_amr_write_2vars.c)
  target_link_libraries(adios_amr_write_2vars adios ${ADIOSLIB_LDADD})
  
  add_executable(adios_staged_read adios_staged_read.c)
  target_link_libraries(adios_staged_read adiosread ${ADIOSREADLIB_LDADD})

  add_executable(adios_staged_read_2vars adios_staged_read_2vars.c)
  target_link_libraries(adios_staged_read_2vars adiosread ${ADIOSREADLIB_LDADD})

  add_executable(write_read write_read.c)
  target_link_libraries(write_read adios ${ADIOSLIB_LDADD})

  add_executable(write_alternate write_alternate.c)
  target_link_libraries(write_alternate adios ${ADIOSLIB_LDADD})

  add_executable(many_vars many_vars.c)
  target_link_libraries(many_vars adios ${ADIOSLIB_LDADD})

  add_executable(selections selections.c)
  target_link_libraries(selections adios ${ADIOSLIB_LDADD})

  add_executable(path_test path_test.c)
  target_link_libraries(path_test adios ${ADIOSLIB_LDADD})

  add_executable(reuse_dim reuse_dim.c)
  target_link_libraries(reuse_dim adios ${ADIOSLIB_LDADD})

  add_executable(adios_transforms_read_write.c)
  target_link_libraries(adios_transforms_read_write adios ${ADIOSLIB_LDADD})
endif(BUILD_WRITE)

if(BUILD_FORTRAN)
  if(BUILD_WRITE)
    add_executable(posix_method posix_method.F90 gwrite_posix_method.fh)
    target_link_libraries(posix_method adiosf ${ADIOSLIB_LDADD})
    include_directories(${PROJECT_BINARY_DIR}/src)
    add_custom_command(
	OUTPUT gwrite_posix_method.fh
	COMMAND ${PROJECT_SOURCE_DIR}/utils/gpp/gpp.py ${PROJECT_SOURCE_DIR}/tests/suite/programs/posix_method.xml
	DEPENDS posix_method.xml
	)
  endif(BUILD_WRITE)
endif(BUILD_FORTRAN)

#install(FILES adios_amr_write_2vars.xml adios_amr_write.xml local_array_time.xml posix_method.xml write_alternate.xml write_read.xml DESTINATION ${PROJECT_BINARY_DIR}/tests/suite/programs)

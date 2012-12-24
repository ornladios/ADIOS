include_directories(${PROJECT_SOURCE_DIR}/tests/bp_read)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/tests/bp_read)

set(genbp_CPPFLAGS "${ADIOSLIB_CPPFLAGS}")
set(genbp_CFLAGS "${ADIOSLIB_CFLAGS}")

set(bp_read_c_CPPFLAGS "-DADIOS_USE_READ_API_1 ${ADIOSLIB_CPPFLAGS}")
set(bp_read_c_CFLAGS "${ADIOSLIB_CFLAGS}")

add_executable(genbp genbp.c)
target_link_libraries(genbp ${ADIOS_LIB} ${ADIOSREAD_LIB} ${ADIOSLIB_LDADD} ${MPI_LIBRARIES})
set_target_properties(genbp PROPERTIES COMPILE_FLAGS "${genbp_CPPFLAGS} ${genbp_CFLAGS}")

add_executable(bp_read_c bp_read_c.c)
target_link_libraries(bp_read_c ${ADIOSREAD_LIB} ${ADIOSREADLIB_LDADD} ${MPI_LIBRARIES})
set_target_properties(bp_read_c PROPERTIES COMPILE_FLAGS "${bp_read_c_CPPFLAGS} ${bp_read_c_CFLAGS}")

if(BUILD_FORTRAN)
  add_executable(bp_read_f bp_read_f.F90)
  target_link_libraries(bp_read_f ${ADIOSREAD_LIB} ${ADIOSREADLIB_LDADD})
  set(bp_read_f_CPPFLAGS "${ADIOSLIB_CPPFLAGS}")
  set(bp_read_f_CFLAGS "${ADIOSLIB_CFLAGS}")
  set_target_properties(bp_read_f PROPERTIES COMPILE_FLAGS "${bp_read_f_CPPFLAGS} ${bp_read_f_CFLAGS}")

  add_executable(bp_read_f_v1 bp_read_f_v1.F90)
  target_link_libraries(bp_read_f_v1 ${ADIOSREAD_V1_LIB} ${ADIOSREADLIB_LDADD})
  set(bp_read_f_v1_CPPFLAGS "${ADIOSLIB_CPPFLAGS}")
  set(bp_read_f_v1_CFLAGS "${ADIOSLIB_CFLAGS}")
  set_target_properties(bp_read_f_v1 PROPERTIES COMPILE_FLAGS "${bp_read_f_v1_CPPFLAGS} ${bp_read_f_v1_CFLAGS}")

endif(BUILD_FORTRAN)



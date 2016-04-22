include_directories(${PROJECT_SOURCE_DIR}/tests/test_src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/public)

include_directories(${PROJECT_BINARY_DIR}/tests/test_src)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)

link_directories(${PROJECT_BINARY_DIR}/tests/test_src)


set(C_PROGS_READONLY hashtest copy_subvolume text_to_pairstruct test_strutil)


if(BUILD_WRITE)
  set(C_PROGS_WRITE transforms_specparse group_free_test)
endif(BUILD_WRITE)


if(BUILD_FORTRAN)
  set(F_PROGS_READONLY selection_api)
  #  if(BUILD_WRITE)
  #    set(F_PROGS_WRITE xxx)
  #  endif(BUILD_WRITE)
endif(BUILD_FORTRAN)

foreach (PROG ${C_PROGS_READONLY} )
  add_executable(${PROG} ${PROG}.c)
  target_link_libraries(${PROG} adiosread_nompi ${ADIOSREADLIB_SEQ_LDADD})
endforeach()

foreach (PROG ${C_PROGS_WRITE} )
  add_executable(${PROG} ${PROG}.c)
  target_link_libraries(${PROG} adios_nompi ${ADIOSLIB_SEQ_LDADD})
endforeach()

foreach (PROG ${F_PROGS_READONLY} )
  add_executable(${PROG} ${PROG}.F90)
  target_link_libraries(${PROG} adiosreadf_nompi ${ADIOSREADLIB_SEQ_LDADD})
endforeach()

foreach (PROG ${F_PROGS_WRITE} )
  add_executable(${PROG} ${PROG}.F90)
  target_link_libraries(${PROG} adiosf_nompi ${ADIOSLIB_SEQ_LDADD})
endforeach()

set_property(TARGET ${C_PROGS_WRITE}
  PROPERTY LINKER_LANGUAGE CXX)
set_property(TARGET ${C_PROGS_READONLY}
  PROPERTY LINKER_LANGUAGE CXX)
set_property(TARGET ${F_PROGS_READONLY}
  PROPERTY LINKER_LANGUAGE CXX)
set_property(TARGET ${F_PROGS_WRITE}
  PROPERTY LINKER_LANGUAGE CXX)

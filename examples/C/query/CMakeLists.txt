include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/query)
link_directories(${PROJECT_BINARY_DIR}/examples/C/query)

add_executable(adios_read_ext_test adios_read_ext_test.c)
target_link_libraries(adios_read_ext_test adios ${ADIOSLIB_LDADD})

if(HAVE_ALACRITY)
  add_executable(adios_build_alac_index adios_build_alac_index.c)
  target_link_libraries(adios_build_alac_index adios ${ADIOSLIB_LDADD})

  add_executable(adios_read_alac_index adios_read_alac_index.c)
  target_link_libraries(adios_read_alac_index adios ${ADIOSLIB_LDADD})

  add_executable(adios_show_bp_file adios_show_bp_file.c)
  target_link_libraries(adios_show_bp_file adios ${ADIOSLIB_LDADD})
endif()

set (PROGS adios_read_ext_test)
if(HAVE_ALACRITY)
  set (PROGS ${PROGS} adios_build_alac_index adios_read_alac_index adios_show_bp_file)
endif()
foreach (PROG ${PROGS} )
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROG} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
endforeach()


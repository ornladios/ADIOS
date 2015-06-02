include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/transforms)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/transforms)

add_executable(transforms_adios_global adios_global.c)
target_link_libraries(transforms_adios_global adios ${ADIOSLIB_LDADD})

add_executable(adios_read_subv adios_read_subv.c)
target_link_libraries(adios_read_subv adiosread ${ADIOSREADLIB_LDADD})

add_executable(adios_read_points adios_read_points.c)
target_link_libraries(adios_read_points adiosread ${ADIOSREADLIB_LDADD})

add_executable(adios_read_wb_subpg adios_read_wb_subpg.c)
target_link_libraries(adios_read_wb_subpg adiosread ${ADIOSREADLIB_LDADD})

set (PROGS transforms_adios_global adios_read_subv adios_read_points adios_read_wb_subpg)
foreach (PROG ${PROGS} )
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROG} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
endforeach()

file(COPY adios_global.xml gread_temperature.ch gwrite_temperature.ch 
     DESTINATION ${PROJECT_BINARY_DIR}/examples/C/transforms)

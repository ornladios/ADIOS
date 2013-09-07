include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/schema)
link_directories(${PROJECT_BINARY_DIR}/examples/C/schema)

add_executable(uniform uniform.c)
target_link_libraries(uniform adios ${ADIOSLIB_LDADD})

add_executable(rectilinear rectilinear.c)
target_link_libraries(rectilinear adios ${ADIOSLIB_LDADD})

add_executable(structured structured.c)
target_link_libraries(structured adios ${ADIOSLIB_LDADD})

add_executable(unstructured unstructured.c)
target_link_libraries(unstructured adios ${ADIOSLIB_LDADD})

set (PROGS uniform rectilinear structured unstructured)
foreach (PROG ${PROGS} )
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROG} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
  target_link_libraries(${PROG} adios ${MPI_C_LIBRARIES})
endforeach()

#install(FILES structured.xml rectilinear.xml unstructured.xml uniform.xml DESTINATION ${PROJECT_BINARY_DIR}/examples/C/schema)

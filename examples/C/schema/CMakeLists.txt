include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/schema)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/schema)

add_executable(uniform2d uniform2d.c)
target_link_libraries(uniform2d adios ${ADIOSLIB_LDADD})

add_executable(rectilinear2d rectilinear2d.c)
target_link_libraries(rectilinear2d adios ${ADIOSLIB_LDADD})

add_executable(structured2d structured2d.c)
target_link_libraries(structured2d adios ${ADIOSLIB_LDADD})

add_executable(tri2d tri2d.c)
target_link_libraries(tri2d adios ${ADIOSLIB_LDADD})

add_executable(triangle2d triangle2d.c)
target_link_libraries(triangle2d adios ${ADIOSLIB_LDADD})

add_executable(uniform2d_noxml uniform2d_noxml.c)
target_link_libraries(uniform2d_noxml adios ${ADIOSLIB_LDADD})

add_executable(rectilinear2d_noxml rectilinear2d_noxml.c)
target_link_libraries(rectilinear2d_noxml adios ${ADIOSLIB_LDADD})

add_executable(structured2d_noxml structured2d_noxml.c)
target_link_libraries(structured2d_noxml adios ${ADIOSLIB_LDADD})

add_executable(tri2d_noxml tri2d_noxml.c)
target_link_libraries(tri2d_noxml adios ${ADIOSLIB_LDADD})

set (PROGS uniform2d rectilinear2d structured2d tri2d triangle2d)
foreach (PROG ${PROGS} )
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROG} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROG} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
  target_link_libraries(${PROG} adios ${MPI_C_LIBRARIES})
endforeach()

file(COPY triangle2d.readme triangle2d.xml triangle2d.n4.png 
        uniform2d.readme uniform2d.xml uniform2d.png
        rectilinear2d.readme rectilinear2d.xml rectilinear2d.png
        structured2d.readme structured2d.xml
        tri2d.readme tri2d.xml tri2d.n.png
        DESTINATION ${PROJECT_BINARY_DIR}/examples/C/schema)

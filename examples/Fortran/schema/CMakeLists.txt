include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/examples/Fortran/schema)
link_directories(${PROJECT_BINARY_DIR}/examples/Fortran/schema)

add_executable(tri2d_f_noxml tri2d_noxml.F90)
target_link_libraries(tri2d_f_noxml adiosf ${ADIOSLIB_LDADD})

add_executable(tri2d_f_noxml_seperate tri2d_noxml_seperate.F90)
target_link_libraries(tri2d_f_noxml_seperate adiosf ${ADIOSLIB_LDADD})

#install(FILES scalars.xml gwrite_scalars.fh gread_scalars.fh job.pbs DESTINATION ${PROJECT_BINARY_DIR}/examples/Fortran/scalars)

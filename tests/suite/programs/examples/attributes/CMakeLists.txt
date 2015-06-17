include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test/suite/programs/examples/attributes)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/test/suite/programs/examples/attributes)

add_executable(attributes_write_C attributes_write_C.c)
if(MPI_COMPILE_FLAGS)
  set_target_properties(attributes_write_C PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(attributes_write_C PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()
target_link_libraries(attributes_write_C adios ${ADIOSLIB_LDADD} ${MPI_C_LIBRARIES})
add_custom_command(
	OUTPUT gwrite_temperature.ch
	COMMAND ${PROJECT_SOURCE_DIR}/utils/gpp/gpp.py ${PROJECT_SOURCE_DIR}/tests/suite/programs/examples/attributes/attributes_C.xml
	DEPENDS attributes_C.xml
	)


add_executable(attributes_read_C attributes_read_C.c)
if(MPI_COMPILE_FLAGS)
  set_target_properties(attributes_read_C PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(attributes_read_C PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()
target_link_libraries(attributes_read_C adiosread ${ADIOSREADLIB_LDADD} ${MPI_C_LIBRARIES})
#set_target_properties(attributes_read PROPERTIES COMPILE_FLAGS "-DADIOS_USE_READ_API_1")

#install(FILES attributes.xml gwrite_temperature.ch gread_temperature.ch DESTINATION ${PROJECT_BINARY_DIR}/examples/C/attributes)


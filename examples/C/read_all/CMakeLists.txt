include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/read_all)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/read_all)

add_executable(read_all read_all.c)
target_link_libraries(read_all adiosread_nompi ${ADIOSLIB_LDADD})
set_target_properties(read_all PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}_NOMPI") 

add_executable(read_all_v1 read_all_v1.c)
target_link_libraries(read_all_v1 adiosread_nompi ${ADIOSREADLIB_LDADD})
set_target_properties(read_all_v1 PROPERTIES COMPILE_FLAGS "${MACRODEFFLAG}_NOMPI ${MACRODEFFLAG}ADIOS_USE_READ_API_1")

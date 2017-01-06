include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/examples/C/icee_arrays/array)
include_directories(${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/examples/C/icee_arrays/array)

add_executable(adios_write adios_write.c)
set_property(TARGET adios_write  PROPERTY COMPILE_FLAGS "${ADIOSLIB_CFLAGS}")
target_link_libraries(adios_write adios ${ADIOSLIB_LDADD})

add_executable(adios_read adios_read.c)
set_property(TARGET adios_read  PROPERTY COMPILE_FLAGS "${ADIOSREADLIB_CFLAGS}")
target_link_libraries(adios_read adios ${ADIOSREADLIB_LDADD})

file(COPY README.txt 
        DESTINATION ${PROJECT_BINARY_DIR}/examples/C/icee_arrays/array)

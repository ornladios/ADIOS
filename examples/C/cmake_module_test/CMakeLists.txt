cmake_minimum_required(VERSION 2.8.0)
PROJECT(schema_uniform)

include_directories(${PROJECT_SOURCE_DIR})
link_directories(${PROJECT_BINARY_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{ADIOS_ROOT}/lib" "$ENV{ADIOS_DIR}/lib" "$ENV{ADIOS_DIR}/etc" "$ENV{INSTALL_PREFIX}/lib")

add_executable(uniform2d uniform2d.c)
find_package(ADIOS)
if(ADIOS_FOUND)
    include_directories(${ADIOS_INCLUDE_DIRS})
    target_link_libraries(uniform2d ${ADIOS_LIBRARIES})
else()
    message("cannot find adios")
endif()

file(COPY uniform2d.xml uniform2d.readme DESTINATION ${PROJECT_BINARY_DIR})

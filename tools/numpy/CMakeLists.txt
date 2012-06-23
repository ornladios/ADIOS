cmake_minimum_required (VERSION 2.8)
set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Modules)
include(${CMAKE_CURRENT_SOURCE_DIR}/Modules/CMakeUtilityFunctions.cmake)

project (AdiosPython CXX)

find_package (PythonLibs REQUIRED 2.7)
include_directories (${PYTHON_INCLUDE_DIRS})
##print(PYTHON_INCLUDE_DIRS)

find_package (PythonInterp REQUIRED)
find_package (Numpy REQUIRED)
include_directories (${PYTHON_NUMPY_INCLUDE_DIR})

option (ADIOS_USE_MPI "Use parallel ADIOS" ON)

find_package (ADIOS REQUIRED)
add_definitions (${ADIOS_COMPILE_FLAGS})  
include_directories (${ADIOS_INCLUDE_PATH})

find_package (MPI REQUIRED)
include_directories (${MPI_CXX_INCLUDE_PATH})

set (PSRCS
  adios.pyx
  )

foreach (FILE ${PSRCS})
  string (REGEX REPLACE ".pyx" ".cpp" OUTFILE ${FILE})

  add_custom_command (
    OUTPUT ${OUTFILE}
    COMMAND cython -o ${CMAKE_CURRENT_BINARY_DIR}/${OUTFILE} -I /ccs/home/jyc/sith/.local/lib/python2.7/site-packages/mpi4py/include ${CMAKE_CURRENT_SOURCE_DIR}/${FILE}
    DEPENDS ${FILE}
    )
  list (APPEND SRCS ${OUTFILE})
endforeach ()

add_library (adios MODULE
  ${SRCS}
  )

target_link_libraries (adios
  ${ADIOS_LIBRARIES}
  ${PYTHON_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  )

set_target_properties (adios
  PROPERTIES OUTPUT_NAME adios
  PREFIX ""
  ##LINK_FLAGS "-undefined dynamic_lookup"
  )

ENABLE_TESTING()
add_test (test_write ${CMAKE_CURRENT_SOURCE_DIR}/test/test_adios_write.sh ${CMAKE_CURRENT_SOURCE_DIR}/test)
add_test (test_read ${CMAKE_CURRENT_SOURCE_DIR}/test/test_adios_read.sh ${CMAKE_CURRENT_SOURCE_DIR}/test)
add_test (test_noxml ${CMAKE_CURRENT_SOURCE_DIR}/test/test_adios_noxml.sh ${CMAKE_CURRENT_SOURCE_DIR}/test)

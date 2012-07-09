cmake_minimum_required (VERSION 2.8)
set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Modules)
include(${CMAKE_CURRENT_SOURCE_DIR}/Modules/CMakeUtilityFunctions.cmake)

project (AdiosJava CXX Java)

include_directories (${CMAKE_CURRENT_BINARY_DIR})

find_package (Java REQUIRED)
include_directories (${JNI_INCLUDE_DIRS})

find_package (JNI REQUIRED)
include_directories (${JNI_INCLUDE_DIRS})

option (ADIOS_USE_MPI "Use parallel ADIOS" ON)

find_package (ADIOS REQUIRED)
add_definitions (${ADIOS_COMPILE_FLAGS})  
add_definitions (-DADIOS_USE_READ_API_1)
include_directories (${ADIOS_INCLUDE_PATH})

if (ADIOS_USE_MPI)
  find_package (MPI REQUIRED)
  include_directories (${MPI_CXX_INCLUDE_PATH})
  add_definitions (-DADIOS_USE_MPI)  
else ()
  find_package (MPI REQUIRED)
endif()

set (SRC
  gov_ornl_ccs_AdiosDatatype.cxx
  gov_ornl_ccs_AdiosFlag.cxx
  gov_ornl_ccs_AdiosBufferAllocWhen.cxx
  gov_ornl_ccs_Adios.cxx
  gov_ornl_ccs_AdiosFile.cxx
  gov_ornl_ccs_AdiosGroup.cxx
  gov_ornl_ccs_AdiosVarinfo.cxx
  )

foreach (FILE ${SRC})
  string (REGEX REPLACE ".cxx" ".h" FILE ${FILE})
  list (APPEND HEADERS ${FILE})
endforeach ()

add_library (AdiosJava SHARED
  ${HEADERS}
  ${SRC}
  )

set (CMAKE_JNI_TARGET TRUE)
set (CMAKE_JAVA_TARGET_VERSION 1.2.0)

target_link_libraries (AdiosJava
  ${ADIOS_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  )

foreach (FILE ${HEADERS})
  string (REGEX REPLACE "_" "/" JAVA_FILE ${FILE})
  string (REGEX REPLACE "\\.h" ".java" JAVA_FILE ${JAVA_FILE})

  string (REGEX REPLACE "_" "." JAVA_OBJ ${FILE})
  string (REGEX REPLACE "\\.h" "" JAVA_OBJ ${JAVA_OBJ})

  add_custom_command (
    OUTPUT ${FILE}
    COMMAND ${Java_JAVAC_EXECUTABLE} -d ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${JAVA_FILE}
    COMMAND ${Java_JAVAH_EXECUTABLE} -jni -d ${CMAKE_CURRENT_BINARY_DIR} ${JAVA_OBJ}
    DEPENDS ${JAVA_FILE} ${PREV}
    )
  list (APPEND PREV ${FILE})
endforeach ()

foreach (FILE ${SRC})
  string (REGEX REPLACE "_" "/" FILE ${FILE})
  string (REGEX REPLACE "\\.cxx" ".java" FILE ${FILE})
  list (APPEND JAVA_SRC ${FILE})
endforeach ()

add_custom_command (TARGET AdiosJava
  POST_BUILD
  COMMAND cmake -E echo "Creating jar file..."
  COMMAND jar cvf AdiosJava.jar gov
  )

install (TARGETS AdiosJava
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/AdiosJava.jar DESTINATION lib)

ENABLE_TESTING()
add_test (test1 ${CMAKE_CURRENT_SOURCE_DIR}/test/test1.sh ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/test AdiosTest)
add_test (test2 ${CMAKE_CURRENT_SOURCE_DIR}/test/test2.sh ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/test AdiosNoxmlTest)


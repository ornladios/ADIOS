find_package(Java REQUIRED)
include(UseJava)

find_package (JNI REQUIRED)
include_directories (${JNI_INCLUDE_DIRS})

option (ADIOS_USE_MPI "Use parallel ADIOS" OFF)

find_package (ADIOS REQUIRED)
add_definitions (${ADIOS_COMPILE_FLAGS})
include_directories (${ADIOS_INCLUDE_PATH})
message (STATUS ADIOS_COMPILE_FLAGS ${ADIOS_COMPILE_FLAGS})

if (ADIOS_USE_MPI)
  find_package (MPI REQUIRED)
  include_directories (${MPI_CXX_INCLUDE_PATH})
  add_definitions (-DADIOS_USE_MPI)
endif()

##find_package (SWIG REQUIRED)
##include(UseSWIG)
##swig_add_module (swig_adios java adios.i)
##swig_link_libraries (swig_adios ${ADIOS_LIBRARIES})

##include_directories (${CMAKE_CURRENT_BINARY_DIR})

if (APPLE)
  set (CMAKE_MACOSX_RPATH ON)
endif (APPLE)

add_library (adiosjava SHARED
  adiosJAVA_wrap.c)

target_link_libraries (adiosjava
  ${ADIOS_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  )

set (CMAKE_JNI_TARGET TRUE)
set (CMAKE_JAVA_TARGET_VERSION 1.8)

file (GLOB SWIG_JAVA_SRC ornl/adios/ext/*.java)
add_jar (ornl-adiosjava
  ornl/adios/Adios.java
  ornl/adios/AdiosFile.java
  ornl/adios/AdiosVar.java
  ornl/adios/AdiosWriter.java
  ${SWIG_JAVA_SRC})

add_custom_target (
    swig
    DEPENDS adios.i
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/swiggen.sh
    COMMENT "Running SWIG to generate files"
)
#add_dependencies (adiosjava genswig)

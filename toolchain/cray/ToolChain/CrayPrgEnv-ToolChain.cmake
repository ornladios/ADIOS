# Make sure we can find the ComputeNodeLinux Platform file since it's not
# shipped with CMake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/..")

set(CMAKE_SYSTEM_NAME ComputeNodeLinux)

# Make sure we have the appropriate environment loaded
set(_CRAYPE_ROOT "$ENV{CRAYPE_DIR}")
if(NOT _CRAYPE_ROOT)
  set(_CRAYPE_ROOT "$ENV{ASYNCPE_DIR}")
endif()
if(NOT _CRAYPE_ROOT)
  message(FATAL_ERROR "Neither the ASYNCPE_DIR or CRAYPE_DIR environment variable are defined but the CrayPrgEnv toolchain module requires one.  This usually means that the necessary PrgEnv-* module is not loaded")
endif()

# Explicitly use the cray compiler wrappers from the PrgEnv-* module
set(CMAKE_C_COMPILER       "${_CRAYPE_ROOT}/bin/cc")
set(CMAKE_CXX_COMPILER     "${_CRAYPE_ROOT}/bin/CC")
set(CMAKE_Fortran_COMPILER "${_CRAYPE_ROOT}/bin/ftn")

# These shouldn't really be necessary since the Cray compiler drivers pay
# attention to the environment variables but this will force the options
# into the generated build files so even if the environment changes after
# configure, the correct commands will still get called.
if(DEFINED ENV{CRAYPE_COMPILE_TARGET})
  list(APPEND _CRAYPE_OPTS "-target=$ENV{CRAYPE_COMPILE_TARGET}")
else()
  if(DEFINED ENV{CRAY_CPU_TARGET})
    list(APPEND _CRAYPE_OPTS "-target-cpu=$ENV{CRAY_CPU_TARGET}")
  endif()
  if(DEFINED ENV{CRAY_ACCEL_TARGET})
    list(APPEND _CRAYPE_OPTS "-target-accel=$ENV{CRAY_ACCEL_TARGET}")
  endif()
  if(DEFINED ENV{CRAY_NETWORK_TARGET})
    list(APPEND _CRAYPE_OPTS "-target-network=$ENV{CRAY_NETWORK_TARGET}")
  endif()
endif()
if("$ENV{CRAYPE_LINK_TYPE}" STREQUAL "dynamic")
  list(APPEND _CRAYPE_OPTS "-dynamic")
else() # Explicit or implicit static
  list(APPEND _CRAYPE_OPTS "-static")
endif()
string(REPLACE ";" " " _CRAYPE_OPTS "${_CRAYPE_OPTS}")

set(CMAKE_C_FLAGS       "${_CRAYPE_OPTS}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS     "${_CRAYPE_OPTS}" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "${_CRAYPE_OPTS}" CACHE STRING "" FORCE)

# Guide the search for binutils. These can't be auto-detected because of
# the forced CMAKE_FIND_ROOT_PATH_MODE_* variables so we explicitly tell
# the CMakeFindBinUtils module to search in the host's directories
set(_CMAKE_TOOLCHAIN_LOCATION /usr/bin)

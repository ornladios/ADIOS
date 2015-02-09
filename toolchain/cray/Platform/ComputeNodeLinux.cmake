# Compute Node Linux doesn't quite work the same as native Linux so all of this
# needs to be custom.  We use the variables defined through Cray's environment
# modules to set up the right paths for things.

set(CMAKE_SYSTEM_VERSION "$ENV{CRAYOS_VERSION}")
if(NOT CMAKE_SYSTEM_VERSION)
  set(CMAKE_SYSTEM_VERSION "$ENV{XTOS_VERSION}")
endif()
if(NOT CMAKE_SYSTEM_VERSION)
  message(FATAL_ERROR "The CMAKE_SYSTEM_VERSION variable is not set and neither the XTOS_VERSION or CRAYOS_VERSION environment variable are defined.  The ComputeNodeLinux CMake platform module either requires it to be manually set or the environment variable to be available. This usually means that the necessary PrgEnv-* module is not loaded")
endif()

# All cray systems are x86 CPUs and have been for quite some time
set(CMAKE_SYSTEM_PROCESSOR "x86_64")

set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
set(CMAKE_STATIC_LIBRARY_PREFIX "lib")
set(CMAKE_STATIC_LIBRARY_SUFFIX ".a")

set(CMAKE_FIND_LIBRARY_PREFIXES "lib")

set(CMAKE_DL_LIBS dl)

# Normally this sort of logic would belong in the toolchain file but the
# order things get loaded in cause anything set here to override the toolchain
# so we'll explicitly check for static compiler options in order to specify
# whether or not the platform will support it.

# If the link type is not explicitly specified in the environment then we'll
# assume that the code will be built statically but that it's dependencies
# can be mixed
if(NOT DEFINED ENV{CRAYPE_LINK_TYPE} OR
   "$ENV{CRAYPE_LINK_TYPE}" STREQUAL "dynamic")
  set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
else() # Explicit static
  set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

# Make sure we have the appropriate environment loaded
if(NOT DEFINED ENV{SYSROOT_DIR})
  #message(WARNING "The SYSROOT_DIR environment variable is not defined.  This is usually due to the appropriate {xe,xt,xc}-sysroot module not being loaded.  A dummy default search path will be used insted, effectively disabling the default search path.")
  set(_CRAY_SYSROOT "/dev/null")
else()
  set(_CRAY_SYSROOT "$ENV{SYSROOT_DIR}")
endif()

include(Platform/UnixPaths)
set(CMAKE_FIND_ROOT_PATH "${_CRAY_SYSROOT}")
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
if(CMAKE_VERSION VERSION_LESS 3.1)
  # This mode is definitely not ideal but necessary for older versions of
  # CMake
  set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY BOTH)
  set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)
  set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE BOTH)
else()
  # Set up system search paths that CMake will use to look for libraries and
  # include files.  These will be the standard UNIX search paths but rooted
  # in the SYSROOT of the compute nodes.  User specified search paths remain
  # untouched
  set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY SYSTEM_ONLY)
  set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE SYSTEM_ONLY)
  set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE SYSTEM_ONLY)
endif()

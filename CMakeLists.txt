cmake_minimum_required(VERSION 2.8)
project(ADIOS_PROJECT)
include_directories(${PROJECT_SOURCE_DIR})
link_directories(${PROJECT_SOURCE_DIR})

#################################### start generating config.h #######################################
INCLUDE (CheckIncludeFiles)
INCLUDE (CheckFunctionExists)

# Define to dummy `main' function (if any) required to link to the Fortran
# libraries.

# Define if F77 and FC dummy `main' functions are identical.

# Define to a macro mangling the given C identifier (in lower and upper case),
# which must not contain underscores, for linking with Fortran.

# As FC_FUNC, but for C identifiers containing underscores.

# Define to 1 if your Fortran compiler doesn't accept -c and -o together.

# Define to 1 if you have the <Attributes.h> header file.
CHECK_INCLUDE_FILES(Attributes.h HAVE_ATTRIBUTES_H)

# Define to 1 if you have the `clock_gettime' function.
CHECK_FUNCTION_EXISTS(clock_gettime HAVE_CLOCK_GETTIME)

# Define if you have CRAY_PMI.
set(HAVE_CRAY_PMI 0)

# Define if you have Cray Portals.
set(HAVE_CRAY_PORTALS 0)

# Define if you have CRAY_UGNI.
set(HAVE_CRAY_UGNI 0)

# Define if you have the DATASPACES.
set(HAVE_DATASPACES 0)

# Define to 1 if you have the <dataspaces.h> header file.
set(HAVE_DATASPACES_H 0)

# Define if you have the DIMES.
set(HAVE_DIMES 0)

# Define to 1 if you have the <dimes.h> header file.
set(HAVE_DIMES_H 0)

# Define to 1 if you have the <dlfcn.h> header file.
CHECK_INCLUDE_FILES(dlfcn.h HAVE_DLFCN_H)

# Define if you have DMALLOC.
set(HAVE_DMALLOC 0)

# Define to 1 if you have the <dmalloc.h> header file.
CHECK_INCLUDE_FILES(dmalloc.h HAVE_DMALLOC_H)

# Define to 1 if you have the `fdatasync' function.
CHECK_FUNCTION_EXISTS(fdatasync HAVE_FDATASYNC)

# Define to 1 if you have the `gettimeofday' function.
CHECK_FUNCTION_EXISTS(gettimeofday HAVE_GETTIMEOFDAY)

# Define if you have HDF5.
set(HAVE_HDF5 0)

# Define to 1 if you have the <hdf5.h> header file.
set(HAVE_HDF5_H 0)

# Define to 1 if you have <infiniband/verbs.h>.
CHECK_INCLUDE_FILES(infiniband/verbs.h HAVE_IBVERBS_H)

# Define if you have the Infiniband.
set(HAVE_INFINIBAND 0)

# Define to 1 if you have the <inttypes.h> header file.
CHECK_INCLUDE_FILES(inttypes.h HAVE_INTTYPES_H)

# Define to 1 if you have the `rt' library (-lrt). 
set(HAVE_LIBRT 0)

# Define if you have LUSTRE.
set(HAVE_LUSTRE 0)

# Define to 1 if you have the <memory.h> header file.
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H)

# Define if you have the MPI library.
set(HAVE_MPI 0)
find_package(MPI)
if(MPI_FOUND)
  set (HAVE_MPI 1)
  message(STATUS "find MPI library")
endif(MPI_FOUND)

# Define if you have the MXML.
# MXML_LIBRARY path /sw/analysis-x64/mxml/2.6/centos5.5_pgi10.9/lib
set(HAVE_MXML 0)
find_library(MXML_LIBRARY NAMES mxml PATHS "")
if(MXML_LIBRARY)
  message(STATUS "find MXML library")
  set(HAVE_MXML 1)
endif(MXML_LIBRARY)

# Define to 1 if you have the <mxml.h> header file.
# mxml.h path /sw/analysis-x64/mxml/2.6/centos5.5_pgi10.9/include
find_path(MXML_INCLUDE_DIR NAMES mxml.h PATHS "")
if(NOT MXML_INCLUDE_DIR_NOTFOUND)
  message(STATUS "find MXML header")
  set(HAVE_MXML_H 1)
  include_directories(${MXML_INCLUDE_DIR})
endif(NOT MXML_INCLUDE_DIR_NOTFOUND)

# Define to 1 if you have the `nanosleep' function.
CHECK_FUNCTION_EXISTS(nanosleep HAVE_NANOSLEEP)

# Define if you have NC4PAR.
set(HAVE_NC4PAR 0)

# Define if you have NETCDF.
set(HAVE_NETCDF 0)

# Define to 1 if you have the <netcdf.h> header file.
set(HAVE_NETCDF_H 0)

# Define if you have NSSI.
set(HAVE_NSSI 0)

# Define to 1 if you have the <p3nal_utcp.h> header file.
set(HAVE_P3NAL_UTCP_H 0)

# Define to 1 if you have the <p3rt/p3rt.h> header file.
set(HAVE_P3RT_P3RT_H 0)

# Define if you have the Pablo.
set(HAVE_PABLO 0)

# Define if you have PHDF5.
# phdf5 library path: /sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4/lib
set(HAVE_PHDF5 0)
find_library(PHDF5_LIBRARY_DIR hdf5 "")
if(PHDF5_LIBRARY_DIR)
  message(STATUS "find PHDF5 library")
  set(HAVE_PHDF5 1)
endif(PHDF5_LIBRARY_DIR)

# Define if you have the Portals.
set(HAVE_PORTALS 0)

# Define to 1 if you have <portals3.h>.
set(HAVE_PORTALS3_H 0)

# Define to 1 if you have <portals/portals3.h>.
set(HAVE_PORTALS_PORTALS3_H 0)

# Define if you have POSIX threads libraries and header files.
set(HAVE_PTHREAD 0)
find_package(Threads)
if(Threads_FOUND)
  set(HAVE_PTHREAD 1)
  #INCLUDE_DIRECTORIES(${Threads_INCLUDE_PATH}) 
endif(Threads_FOUND)

# Define to 1 if you have the `pthread_yield' function.
CHECK_FUNCTION_EXISTS(pthread_yield HAVE_PTHREAD_YIELD)

# Define to 1 if you have the `PtlACEntry' function.
CHECK_FUNCTION_EXISTS(PtlACEntry HAVE_PTLACENTRY)
#set(HAVE_PTLACENTRY 0)

# Define to 1 if you have the `PtlErrorStr' function.
CHECK_FUNCTION_EXISTS(PtlErrorStr HAVE_PTLERRORSTR)

# Define to 1 if you have the `PtlEventKindStr' function.
CHECK_FUNCTION_EXISTS(PtlEventKindStr HAVE_PTLEVENTKINDSTR)

# Define to 1 if you have the `PtlGetJid' function.
CHECK_FUNCTION_EXISTS(PtlGetJid HAVE_PTLGETJID)

# Define to 1 if you have the `PtlNIFailStr' function.
CHECK_FUNCTION_EXISTS(PtlNIFailStr HAVE_PTLNIFAILSTR)

# Define to 1 if the system has the type `ptl_eq_handler_t'.
set(HAVE_PTL_EQ_HANDLER_T 0)

# Define if you have PTL_NOACK_REQ.
set(HAVE_PTL_NOACK_REQ 0)

# Define if you have PTL_NO_ACK_REQ.
set(HAVE_PTL_NO_ACK_REQ 0)

# Define to 1 if the system has the type `ptl_time_t'.
set(HAVE_PTL_TIME_T 0)

# Define to 1 if you have the `sched_yield' function.
CHECK_FUNCTION_EXISTS(sched_yield HAVE_SCHED_YIELD)

# Define to 1 if you have the <stdint.h> header file.
CHECK_INCLUDE_FILES(stdint.h HAVE_STDINT_H)

# Define to 1 if you have the <stdlib.h> header file.
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)

# Define to 1 if you have the `strerror' function.
CHECK_FUNCTION_EXISTS(strerror HAVE_STRERROR)

# Define to 1 if cpp supports the ANSI  stringizing operator.
set(HAVE_STRINGIZE 0)

# Define to 1 if you have the <strings.h> header file.
CHECK_INCLUDE_FILES(strings.h HAVE_STRINGS_H)

# Define to 1 if you have the <string.h> header file.
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)

# Define to 1 if you have the `strncpy' function.
CHECK_FUNCTION_EXISTS(strncpy HAVE_STRNCPY)

# Define to 1 if you have the <sys/stat.h> header file.
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H)

# Define to 1 if you have the <sys/types.h> header file.
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)

# Define to 1 if you have the <unistd.h> header file.
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)

# Datatap is disabled
set(NO_DATATAP 1)

# Define to 1 if your C compiler doesn't accept -c and -o together.

# Name of package
set(PACKAGE adios)

# Define to the address where bug reports for this package should be sent.

# Define to the full name of this package.

# Define to the full name and version of this package.

# Define to the one symbol short name of this package.

# Define to the version of this package.

# Define to the necessary symbol if this constant uses a non-standard name on
# your system.

# Skel timing is enable

# Define to 1 if you have the ANSI C header files.

# Need to use MPI compiler for sequential utils

# Major version number
set(VERSION_MAJOR 1)

# Minor version number
set(VERSION_MINOR 4)

# Micro version number
set(VERSION_MICRO 1)

# Version number of package
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_MICRO})

# svn revision number
set(VERSION_SVNREV 2007)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${PROJECT_SOURCE_DIR}/config.h)
#################################### end of generating config.h #######################################

add_definitions(-DO_LARGEFILE=0)

add_subdirectory(src)
#add_subdirectory(utils)
#add_subdirectory(tests)
#add_subdirectory(examples)

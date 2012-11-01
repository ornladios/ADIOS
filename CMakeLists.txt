cmake_minimum_required(VERSION 2.8)
project(ADIOS_PROJECT)
#set(LANGUAGE "" CACHE FILEPATH "")  #C/CXX/Fortran
#enable_language(${LANGUAGE})
#message(STATUS "CMAKE_${LANGUAGE}_COMPILER_WORKS = ${CMAKE_${LANGUAGE}_COMPILER_WORKS}")
include_directories(${PROJECT_SOURCE_DIR})
link_directories(${PROJECT_BINARY_DIR})

set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${prefix})
set(bindir ${exec_prefix}/bin)
set(datadir ${prefix}/share)
set(libdir ${exec_prefix}/lib)
set(includedir ${prefix}/include)

option(BUILD_WRITE "BUILD_WRITE" OFF)
option(BUILD_FORTRAN "BUILD_FORTRAN" OFF)

set(LIBS "")
set(CPPFLAGS "-DMPICH_IGNORE_CXX_SEEK")
set(CFLAGS "-g -fPIC")

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
  if(MPI_LIBRARIES AND MPI_INCLUDE_PATH)
    message(STATUS "find MPI library and header")
    set (HAVE_MPI 1)
    set(MPILIBS ${MPI_LIBRARIES})
    include_directories(${MPI_INCLUDE_PATH})
    set(LIBS "${LIBS} ${MPILIBS}")  
endif(MPI_LIBRARIES AND MPI_INCLUDE_PATH)
endif(MPI_FOUND)

# Define to 1 if you have the <mxml.h> header file.
find_path(MXML_INCLUDE_DIR NAMES mxml.h PATHS ${MXML_DIR}/include)
if(NOT MXML_INCLUDE_DIR_NOTFOUND)
  message(STATUS "find MXML header")
  set(HAVE_MXML_H 1)
  include_directories(${MXML_INCLUDE_DIR})
  set(MXML_CPPFLAGS ${MXML_INCLUDE_DIR})
endif(NOT MXML_INCLUDE_DIR_NOTFOUND)

# Define if you have the MXML.
# MXML: /sw/analysis-x64/mxml/2.6/centos5.5_pgi10.9
set(HAVE_MXML 0)
set(MXML_DIR "" CACHE FILEPATH "") 
find_library(MXML_LIB_DIR NAMES mxml PATHS ${MXML_DIR}/lib)
if(MXML_LIB_DIR AND HAVE_MXML_H)
  message(STATUS "find MXML library")
  set(HAVE_MXML 1)
  set(MXML_LIB ${MXML_LIB_DIR})
  set(LIBS "${LIBS} ${MXML_LIB}")
  set(CPPFLAGS "${CPPFLAGS} ${MXML_CPPFLAGS}")
endif(MXML_LIB_DIR AND HAVE_MXML_H)

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
# phdf5: /sw/sith/hdf5/1.8.5/centos5.5_pgi10.9_ompi1.4.4
set(HAVE_PHDF5 0)
set(PHDF5_DIR "" CACHE FILEPATH "")
find_library(PHDF5_LIB_DIR NAMES hdf5 PATHS "${PHDF5_DIR}/lib")
if(PHDF5_LIB_DIR)
  message(STATUS "find PHDF5 library")
  set(HAVE_PHDF5 1)
endif(PHDF5_LIB_DIR)

# Define if you have the Portals.
set(HAVE_PORTALS 0)

# Define to 1 if you have <portals3.h>.
#set(HAVE_PORTALS3_H 0)
CHECK_INCLUDE_FILES(portals3.h HAVE_PORTALS3_H)

# Define to 1 if you have <portals/portals3.h>.
#set(HAVE_PORTALS_PORTALS3_H 0)
CHECK_INCLUDE_FILES(portals/portals3.h HAVE_PORTALS_PORTALS3_H)

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
################################ end of generating config.h ###############################

set(ADIOSLIB_LDADD "m ${MXML_LIB}")
set(ADIOSLIB_CFLAGS "")
set(ADIOSLIB_CPPFLAGS ${MXML_CPPFLAGS})
set(ADIOSLIB_EXTRA_CPPFLAGS "-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_GNU_SOURCE")
set(ADIOSLIB_SEQ_CPPFLAGS "-D_NOMPI ${MXML_CPPFLAGS}")
set(ADIOSLIB_SEQ_CFLAGS "")
set(ADIOSLIB_SEQ_LDADD "m ${MXML_LIB}")
set(ADIOSLIB_INT_CPPFLAGS "-D_NOMPI ${MXML_CPPFLAGS}")
set(ADIOSLIB_INT_CFLAGS "")
set(ADIOSLIB_INT_LDADD "m ${MXML_LIB}")
set(ADIOSREADLIB_CPPFLAGS "")
set(ADIOSREADLIB_CFLAGS "")
set(ADIOSREADLIB_LDADD "m")
set(ADIOSREADLIB_SEQ_CPPFLAGS "-D_NOMPI")
set(ADIOSREADLIB_SEQ_CFLAGS "")
set(ADIOSREADLIB_SEQ_LDADD "m")

set(NSSI_SERVER_LIBS "")
############################### srart of top CMakeLists.txt ###############################
install(PROGRAMS adios_config DESTINATION ${bindir})

if(HAVE_NSSI)
install(PROGRAMS scripts/create.nssi.config.sh DESTINATION ${bindir}) 
install(PROGRAMS scripts/kill.nssi.staging.sh DESTINATION ${bindir})
install(PROGRAMS scripts/start.nssi.staging.sh DESTINATION  ${bindir})
endif(HAVE_NSSI)

add_definitions(-DO_LARGEFILE=0)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(utils)
add_subdirectory(examples)
############################### end of top CMakeLists.txt ###############################

cmake_minimum_required(VERSION 2.8.6)
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

#option(PHDF5 "PHDF5" OFF)
#option(LUSTRE "LUSTRE" OFF)
if(MACRODEF)
  set(MACRODEFFLAG "-WL,-D")
else()
  set(MACRODEFFLAG "-D")
endif()
message(STATUS "MECRODEFFLAG = ${MACRODEFFLAG}")

set(LIBS "")
set(CPPFLAGS $ENV{CPPFLAGS})
set(CFLAGS $ENV{CFLAGS})
add_definitions(-DO_LARGEFILE=0)

#################################### start generating config.h #######################################
set(BUILD_WITH_CMAKE 1)

INCLUDE (CheckIncludeFiles)
INCLUDE (CheckFunctionExists)

# Define to dummy `main' function (if any) required to link to the Fortran
# libraries.

# Define if F77 and FC dummy `main' functions are identical.

# Define to a macro mangling the given C identifier (in lower and upper case),
# which must not contain underscores, for linking with Fortran.
set(FC_FUNC FC_GLOBAL)

# As FC_FUNC, but for C identifiers containing underscores.
set(FC_FUNC_ FC_GLOBAL_)
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

# Define if you have the MPI library.
if(BUILD_FORTRAN)
  find_package(MPI COMPONENTS Fortran)
  if(MPI_FOUND)
    if(MPI_LIBRARIES AND MPI_INCLUDE_PATH)
      message(STATUS "find MPI library and header")
      set (HAVE_MPI 1)
      set(MPILIBS ${MPI_LIBRARIES})
      include_directories(${MPI_INCLUDE_PATH})
      set(LIBS ${MPILIBS})
    else()
      if(CMAKE_Fortran_COMPILER)
      set (HAVE_MPI 1)
      endif()
    endif(MPI_LIBRARIES AND MPI_INCLUDE_PATH)
  else()
      if(CMAKE_Fortran_COMPILER)
      set (HAVE_MPI 1)
      endif()
  endif(MPI_FOUND)
else()
  find_package(MPI)
  if(MPI_FOUND)
#    message(STATUS "${MPI_LIBRARIES}")
    if(MPI_LIBRARIES AND MPI_INCLUDE_PATH) 
      message(STATUS "find MPI library and header")
      set (HAVE_MPI 1)
      set(MPILIBS ${MPI_LIBRARIES})
#    set(MPILIBS "-lmpi")
      include_directories(${MPI_INCLUDE_PATH})
      set(LIBS ${MPILIBS})
    else()
      if(CMAKE_CXX_COMPILER OR CMAKE_CXX_COMPILER)
      set (HAVE_MPI 1)
      endif()
    endif(MPI_LIBRARIES AND MPI_INCLUDE_PATH)
  else()
    if(CMAKE_CXX_COMPILER OR CMAKE_CXX_COMPILER)
      set (HAVE_MPI 1)
    endif() 
  endif(MPI_FOUND)
endif(BUILD_FORTRAN)

# Define to 1 if you have the <hdf5.h> header file.
if(HDF5)
# Define if you have HDF5.
  if(HDF5_LIBS AND HDF5_FLAGS)
    set(HAVE_HDF5 1)
    find_path(HDF5_INCLUDE_DIR hdf5.h ${HDF5_DIR}/include)
    if(HDF5_INCLUDE_DIR)
      set(HAVE_HDF5_H 1)
      include_directories(${HDF5_INCLUDE_DIR})
      set(CPPFLAGS "${CPPFLAGS} ${HDF5_FLAGS}")
      set(LIBS "${LIBS} ${HDF5_LIBS}")
    else()
      set(HAVE_HDF5_H 0)
    endif()
  else()
    set(HAVE_HDF5 0)
  endif()
else()
  set(HAVE_HDF5 0)
  set(HAVE_HDF5_H 0)
endif(HDF5)

# Define to 1 if you have <infiniband/verbs.h>.
CHECK_INCLUDE_FILES(infiniband/verbs.h HAVE_IBVERBS_H)

# Define if you have the Infiniband.
set(HAVE_INFINIBAND 0)

# Define to 1 if you have the <inttypes.h> header file.
CHECK_INCLUDE_FILES(inttypes.h HAVE_INTTYPES_H)

# Define to 1 if you have the `rt' library (-lrt). 
set(HAVE_LIBRT 0)

# Define if you have LUSTRE.
if(LUSTRE)
#  set(LUSTRE_DIR "" CACHE FILEPATH "")
  if(LUSTRE_DIR)
    find_library(LUSTRE_LIB_DIR NAMES lustre PATHS "${LUSTRE_DIR}/lib")
    if(LUSTRE_LIB_DIR AND ${LUSTRE_DIR}/include)
      message(STATUS "find LUSTRE library")
      set(HAVE_LUSTRE 1)
      set(LIBS "${LIBS} -llustreapi")
      set(LUSTRE_CPPFLAGS "-I${LUSTRE_DIR}/include")
      set(CPPFLAGS "${CPPFLAGS} ${LUSTRE_CPPFLAGS}")
    endif()
  else()
    set(HAVE_LUSTRE 0)
  endif(LUSTRE_DIR)
else()
  set(HAVE_LUSTRE 0)
endif(LUSTRE)

# Define to 1 if you have the <memory.h> header file.
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H)

# Define to 1 if you have the <mxml.h> header file.
find_path(MXML_INCLUDE_DIR NAMES mxml.h PATHS ${MXML_DIR}/include)
if(NOT MXML_INCLUDE_DIR-NOTFOUND)
  message(STATUS "find mxml header")
  set(HAVE_MXML_H 1)
  include_directories(${MXML_INCLUDE_DIR})
  set(MXML_CPPFLAGS "-I${MXML_INCLUDE_DIR}")
  set(CPPFLAGS "${CPPFLAGS} ${MXML_CPPFLAGS}")
endif(NOT MXML_INCLUDE_DIR-NOTFOUND)

# Define if you have the MXML.
# MXML: /sw/analysis-x64/mxml/2.6/centos5.5_pgi10./9
set(HAVE_MXML 0)
find_library(MXML_LIBS NAMES mxml PATHS $ENV{MXML_DIR}/lib)
if(MXML_LIBS)
  message(STATUS "find mxml library")
endif()
if(MXML_LIBS AND HAVE_MXML_H)
  set(HAVE_MXML 1)
  set(LIBS "${LIBS} -lmxml")
  set(CPPFLAGS "${CPPFLAGS} ${MXML_CPPFLAGS}")
endif()

# Define to 1 if you have the `nanosleep' function.
CHECK_FUNCTION_EXISTS(nanosleep HAVE_NANOSLEEP)

# Define if you have NETCDF.
set(HAVE_NETCDF 0)

# Define to 1 if you have the <netcdf.h> header file.
set(HAVE_NETCDF_H 0)

# Define if you have NC4PAR.
set(HAVE_NC4PAR 0)

# Use NC4PAR for bp2ncd if no sequential NETCDF was found
if(HAVE_NC4PAR)
  if(NOT HAVE_NETCDF)
    set(NETCDF_CPPFLAGS "${NC4PAR_CPPFLAGS}")
    set(NETCDF_CFLAGS "${NC4PAR_CFLAGS}")
    set(NETCDF_LIBS "${NC4PAR_LIBS}")
    set(HAVE_NETCDF 1)
    set(USE_PARALLEL_COMPILER 1)
    message(STATUS "Use NC4PAR for sequential NETCDF tools")
  endif()
endif()

# Define if you have NSSI.
set(HAVE_NSSI 0)

# Define to 1 if you have the <p3nal_utcp.h> header file.
set(HAVE_P3NAL_UTCP_H 0)

# Define to 1 if you have the <p3rt/p3rt.h> header file.
set(HAVE_P3RT_P3RT_H 0)

# Define if you have the Pablo.
set(HAVE_PABLO 0)

# Define if you have PHDF5.
if(PHDF5)
  if(PHDF5_LIBS AND PHDF5_FLAGS)
    set(HAVE_PHDF5 1)  
#  find_library(PHDF5_LIB_DIR NAMES hdf5 hdf5_hl PATHS "${PHDF5_DIR}/lib")
    find_path(PHDF5_INCLUDE_DIR hdf5.h ${PHDF5_DIR}/include)
    if(PHDF5_INCLUDE_DIR)
      set(HAVE_PHDF5_H 1)
      include_directories(${PHDF5_INCLUDE_DIR})
      set(CPPFLAGS "${CPPFLAGS} ${PHDF5_FLAGS}")
      set(LIBS "${LIBS} ${PHDF5_LIBS}")
    endif()
  else()
    set(HAVE_PHDF5 0)
  endif()
else()
  set(HAVE_PHDF5 0)
endif(PHDF5)

# Use PHDF5 for bp2h5 if no sequential HDF5 was found
if(HAVE_PHDF5)
  if(NOT HAVE_HDF5)
    set(HDF5_FLAGS "${PHDF5_FLAGS}")
    set(HDF5_LIBS "${PHDF5_LIBS}")
    set(HAVE_HDF5 1)
#    include_directories(${PHDF5_INCLUDE_DIR})
    set(USE_PARALLEL_COMPILER 1)
    set(CPPFLAGS "${CPPFLAGS} ${PHDF5_FLAGS}")
    set(LIBS "${LIBS} ${PHDF5_LIBS}")
    message(STATUS "Use PHDF5 for sequential HDF5 tools")
  endif()
endif()


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

################################ start of Fortran processing ##############################
if(BUILD_FORTRAN)
enable_language(Fortran)
message(STATUS "CMAKE_Fortran_COMPILER_WORKS =${CMAKE_Fortran_COMPILER_WORKS}")

include(FortranCInterface)
FortranCInterface_HEADER(${CMAKE_CURRENT_SOURCE_DIR}/FC.h MACRO_NAMESPACE "FC_")
#FortranCInterface_VERIFY(CXX)
endif()
################################ end of Fortran processing ##############################


# Generate flags for adios lib creation which depends on the packages we found
# ADIOSLIB_* for compiling libadios.a and libadiosf.a, and linking codes with those
# ADIOSLIB_SEQ_* for compiling libadios_nompi.a and libadiosf_nompi.a with dummy MPI
# ADIOSLIB_INT_* for codes using libadios_internal_nompi.a (like bp2ncd, bp2h5, bpdump etc)
# ADIOSREADLIB_* for codes using libadiosread.a, libadiosreadf.a
# ADIOSREADLIB_SEQ_* for codes using libadiosread_nompi.a
# adios_config will contain flags for users' code linking based on these flags
# EXTRA flags are used for compiling adios stuff, which need not go into adios_config for user codes
set(ADIOSLIB_EXTRA_CPPFLAGS "${MACRODEFFLAG}_FILE_OFFSET_BITS=64 ${MACRODEFFLAG}_LARGEFILE_SOURCE ${MACRODEFFLAG}_LARGEFILE64_SOURCE ${MACRODEFFLAG}_GNU_SOURCE")
set(ADIOSLIB_LDADD "m ${MXML_LIBS}")
set(ADIOSLIB_CFLAGS "")
set(ADIOSLIB_CPPFLAGS "${MXML_CPPFLAGS}")
set(ADIOSLIB_SEQ_CPPFLAGS "${MACRODEFFLAG}_NOMPI ${MXML_CPPFLAGS}")
set(ADIOSLIB_SEQ_CFLAGS "")
set(ADIOSLIB_SEQ_LDADD "m ${MXML_LIBS}")
set(ADIOSLIB_INT_CPPFLAGS "${MACRODEFFLAG}_NOMPI ${MXML_CPPFLAGS}")
set(ADIOSLIB_INT_CFLAGS "")
set(ADIOSLIB_INT_LDADD "m ${MXML_LIBS}")
set(ADIOSREADLIB_CPPFLAGS "")
set(ADIOSREADLIB_CFLAGS "")
set(ADIOSREADLIB_LDADD "m")
set(ADIOSREADLIB_SEQ_CPPFLAGS "${MACRODEFFLAG}_NOMPI")
set(ADIOSREADLIB_SEQ_CFLAGS "")
set(ADIOSREADLIB_SEQ_LDADD "m")

if(NOT NO_DATATAP)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${DT_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${DT_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${DT_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${DT_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${DT_CFLAGS}")
  set(ADIOSREADLIB_LDADD "${ADIOSREADLIB_LDADD} ${DT_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${DT_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${DT_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${DT_LIBS}")
endif()

if(HAVE_DATASPACES)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${DATASPACES_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${DATASPACES_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${DATASPACES_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${DATASPACES_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${DATASPACES_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${DATASPACES_CFLAGS}")
  set(ADIOSREADLIB_LDADD "${ADIOSREADLIB_LDADD} ${DATASPACES_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${DATASPACES_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${DATASPACES_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${DATASPACES_LIBS}")
endif()

if(HAVE_DIMES)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${DIMES_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${DIMES_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${DIMES_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${DIMES_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${DIMES_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${DIMES_CFLAGS}")
  set(ADIOSREADLIB_LDADD "${ADIOSREADLIB_LDADD} ${DIMES_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${DIMES_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${DIMES_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${DIMES_LIBS}")
endif()

if(HAVE_NSSI)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${NSSI_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${NSSI_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${NSSI_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${NSSI_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS} ${NSSI_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${NSSI_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${MXML_CPPFLAGS} ${NSSI_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${MXML_CPPFLAGS} ${NSSI_CFLAGS}")
  set(ADIOSREADLIB_LDADD "${ADIOSREADLIB_LDADD} ${MXML_LIBS} ${NSSI_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${MXML_CPPFLAGS} ${NSSI_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${MXML_CPPFLAGS} ${NSSI_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${MXML_LIBS} ${NSSI_LIBS}")
endif()

if(HAVE_NC4PAR)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${NC4PAR_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${NC4PAR_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${NC4PAR_LIBS}")
endif()

if(HAVE_PHDF5)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${PHDF5_FLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${PHDF5_FLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${PHDF5_LIBS}")
  message(STATUS "Use PHDF5 for sequential HDF5 tools")
endif()

if(HAVE_PORTALS)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${PORTALS_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${PORTALS_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${PORTALS_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${PORTALS_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS} ${PORTALS_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${PORTALS_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${PORTALS_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${PORTALS_CFLAGS}")
  set(ADIOSREADLIB_LDADD "${ADIOSREADLIB_LDADD} ${PORTALS_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${PORTALS_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${PORTALS_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${PORTALS_LIBS}")
endif()

if(HAVE_INFINIBAND)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${INFINIBAND_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${INFINIBAND_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${INFINIBAND_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${INFINIBAND_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS} ${INFINIBAND_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${INFINIBAND_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${INFINIBAND_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${INFINIBAND_CFLAGS}")
  set(ADIOSREADLIB_LDADD "${ADIOSREADLIB_LDADD} ${INFINIBAND_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${INFINIBAND_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${INFINIBAND_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${INFINIBAND_LIBS}")
endif()

if(HAVE_CRAY_PMI)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${CRAY_PMI_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${CRAY_PMI_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${CRAY_PMI_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${CRAY_PMI_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS} ${CRAY_PMI_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${CRAY_PMI_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${CRAY_PMI_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${CRAY_PMI_CFLAGS}")
  set(ADIOSREADLIB_LDADD "${ADIOSREADLIB_LDADD} ${CRAY_PMI_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${CRAY_PMI_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${CRAY_PMI_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${CRAY_PMI_LIBS}")
endif()

if(HAVE_CRAY_UGNI)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${CRAY_UGNI_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${CRAY_UGNI_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${CRAY_UGNI_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${CRAY_UGNI_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS} ${CRAY_UGNI_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${CRAY_UGNI_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${CRAY_UGNI_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${CRAY_UGNI_CFLAGS}")
  set(ADIOSREADLIB_LDADD="${ADIOSREADLIB_LDADD} ${CRAY_UGNI_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${CRAY_UGNI_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${CRAY_UGNI_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${CRAY_UGNI_LIBS}")
endif()

if(HAVE_NETCDF)
  set(ADIOSLIB_INT_CPPFLAGS "${ADIOSLIB_INT_CPPFLAGS} ${NETCDF_CPPFLAGS}")
  set(ADIOSLIB_INT_CFLAGS "${ADIOSLIB_INT_CFLAGS} ${NETCDF_CFLAGS}")
  set(ADIOSLIB_INT_LDADD "${ADIOSLIB_INT_LDADD} ${NETCDF_LIBS}")
endif()

if(HAVE_HDF5)
  set(ADIOSLIB_INT_CPPFLAGS "${ADIOSLIB_INT_CPPFLAGS} ${HDF5_CPPFLAGS}")
  set(ADIOSLIB_INT_CFLAGS "${ADIOSLIB_INT_CFLAGS} ${HDF5_CFLAGS}")
  set(ADIOSLIB_INT_LDADD "${ADIOSLIB_INT_LDADD} ${HDF5_LIBS}")
endif()

if(HAVE_DMALLOC)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${MACRODEFFLAG}DMALLOC ${DMALLOC_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${DMALLOC_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${DMALLOC_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${DMALLOC_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS} ${DMALLOC_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${DMALLOC_LIBS}")
  set(ADIOSLIB_INT_CPPFLAGS "${ADIOSLIB_INT_CPPFLAGS} ${DMALLOC_CPPFLAGS}")
  set(ADIOSLIB_INT_CFLAGS "${ADIOSLIB_INT_CFLAGS} ${DMALLOC_CFLAGS}")
  set(ADIOSLIB_INT_LDADD "${ADIOSLIB_INT_LDADD} ${DMALLOC_LIBS}")
  set(ADIOSREADLIB_CPPFLAGS "${ADIOSREADLIB_CPPFLAGS} ${DMALLOC_CPPFLAGS}")
  set(ADIOSREADLIB_CFLAGS "${ADIOSREADLIB_CFLAGS} ${DMALLOC_CFLAGS}")
  set(ADIOSREADLIB_LDADD "${ADIOSREADLIB_LDADD} ${DMALLOC_LIBS}")
  set(ADIOSREADLIB_SEQ_CPPFLAGS "${ADIOSREADLIB_SEQ_CPPFLAGS} ${DMALLOC_CPPFLAGS}")
  set(ADIOSREADLIB_SEQ_CFLAGS "${ADIOSREADLIB_SEQ_CFLAGS} ${DMALLOC_CFLAGS}")
  set(ADIOSREADLIB_SEQ_LDADD "${ADIOSREADLIB_SEQ_LDADD} ${DMALLOC_LIBS}")
endif()

if(HAVE_LUSTRE)
  set(ADIOSLIB_CPPFLAGS "${ADIOSLIB_CPPFLAGS} ${LUSTRE_CPPFLAGS}")
  set(ADIOSLIB_CFLAGS "${ADIOSLIB_CFLAGS} ${LUSTRE_CFLAGS}")
  set(ADIOSLIB_LDADD "${ADIOSLIB_LDADD} ${LUSTRE_LIBS}")
  set(ADIOSLIB_SEQ_CPPFLAGS "${ADIOSLIB_SEQ_CPPFLAGS} ${LUSTRE_CPPFLAGS}")
  set(ADIOSLIB_SEQ_CFLAGS "${ADIOSLIB_SEQ_CFLAGS} ${LUSTRE_CFLAGS}")
  set(ADIOSLIB_SEQ_LDADD "${ADIOSLIB_SEQ_LDADD} ${LUSTRE_LIBS}")
endif()

#set(NSSI_SERVER_LIBS "")
############################### srart of top CMakeLists.txt ###############################
#install(PROGRAMS adios_config DESTINATION ${bindir})

#if(HAVE_NSSI)
#install(PROGRAMS scripts/create.nssi.config.sh DESTINATION ${bindir}) 
#install(PROGRAMS scripts/kill.nssi.staging.sh DESTINATION ${bindir})
#install(PROGRAMS scripts/start.nssi.staging.sh DESTINATION  ${bindir})
#endif(HAVE_NSSI)



#add_definitions(-DO_LARGEFILE=0)

add_subdirectory(src)
#set(ADIOS_LIB ${PROJECT_BINARY_DIR}/src/libadios.a)
#set(ADIOS_LIB adios)
#set(ADIOS_INTERNAL_NOMPI_LIB ${PROJECT_BINARY_DIR}/src/libadios_internal_nompi.a)
#set(ADIOS_INTERNAL_NOMPI_LIB adios_internal_nompi)
#set(ADIOSREAD_NOMPI_LIB ${PROJECT_BINARY_DIR}/src/libadiosread_nompi.a)
#set(ADIOSREAD_NOMPI_LIB adiosread_nompi)
#if(BUILD_FORTRAN)
#  set(ADIOSREAD_LIB ${PROJECT_BINARY_DIR}/src/libadiosreadf.a)
#  set(ADIOSREAD_LIB adiosreadf)
#  set(ADIOSREAD_V1_LIB ${PROJECT_BINARY_DIR}/src/libadiosreadf_v1.a)
#  set(ADIOSREAD_V1_LIB adiosreadf_v1)
#else()
#  set(ADIOSREAD_LIB ${PROJECT_BINARY_DIR}/src/libadiosread.a)
#  set(ADIOSREAD_LIB adiosread)
#endif()

get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
message(STATUS "${_LANGUAGES_}")

add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(examples)
#set_source_files_properties
############################### end of top CMakeLists.txt ###############################

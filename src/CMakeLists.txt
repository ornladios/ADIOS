include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)

if(HAVE_NSSI)
  add_executable(nssi-staging-server nssi/nssi_staging_server.cpp nssi/aggregation.cpp)
  set_target_properties(nssi-staging-server PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS} -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL") 
  target_link_libraries(nssi-staging-server ${PROJECT_BINARY_DIR}/src/libadios.a ${ADIOSLIB_LDADD} ${NSSI_SERVER_LIBS})

  add_executable(nssi-coupling-server nssi/nssi_coupling_server.cpp)
  set_target_properties(nssi-coupling-server PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS} -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL") 
  target_link_libraries(nssi-coupling-server ${PROJECT_BINARY_DIR}/src/libadios.a ${ADIOSLIB_LDADD} ${NSSI_SERVER_LIBS})
endif(HAVE_NSSI)

set(include_HEADERS "")
set(nodist_include_HEADERS "")
set(lib_LIBRARIES "")

set(libadiosread_a_CPPFLAGS "")
set(libadiosreadf_a_CPPFLAGS "")
set(libadiosread_nompi_a_CPPFLAGS "")
set(libadiosreadf_nompi_a_CPPFLAGS "")

set(libadiosread_a_CFLAGS "")
set(libadiosreadf_a_CFLAGS "")
set(libadiosread_nompi_a_CFLAGS "")
set(libadiosreadf_nompi_a_CFLAGS "")

if(BUILD_WRITE)
  set(libadios_a_SOURCES core/adios.c 
                     core/common_adios.c
                     core/adios_internals.c
                     core/adios_internals_mxml.c 
                     core/buffer.c 
                     core/adios_bp_v1.c  
                     core/adios_endianness.c 
                     core/bp_utils.c 
                     core/futils.c 
                     core/adios_error.c 
                     core/adios_read.c 
                     core/adios_read_v1.c 
                     core/common_read.c 
                     core/globals.c 
                     core/adios_timing.c 
                     core/adios_read_hooks.c 
                     core/adios_transport_hooks.c 
                     core/adios_socket.c 
                     core/adios_logger.c 
                     core/util.c 
                     read/read_bp.c 
                     read/read_bp_staged.c 
                     read/read_bp_staged1.c
                     write/adios_mpi.c
                     write/adios_mpi_lustre.c
                     write/adios_mpi_amr.c
                     write/adios_posix.c
                     write/adios_posix1.c)

set(libadios_nompi_a_SOURCES core/adios.c 
                     core/common_adios.c 
                     core/adios_internals.c 
                     core/adios_internals_mxml.c 
                     core/buffer.c 
                     core/adios_bp_v1.c  
                     core/adios_endianness.c 
                     core/bp_utils.c 
                     core/futils.c 
                     core/adios_error.c 
                     core/adios_read.c 
                     core/adios_read_v1.c 
                     core/common_read.c 
                     core/globals.c 
                     core/mpidummy.c 
                     core/adios_timing.c 
                     core/adios_read_hooks.c 
                     core/adios_transport_hooks.c 
                     core/adios_socket.c 
                     core/adios_logger.c 
                     core/util.c 
                     read/read_bp.c 
                     read/read_bp_staged.c 
                     read/read_bp_staged1.c 
                     write/adios_posix.c 
                     write/adios_posix1.c)

  if(BUILD_FORTRAN)
    set(FortranLibSources core/adiosf.c 
                       core/common_adios.c 
                       core/adios_internals.c 
                       core/adios_internals_mxml.c
                       core/buffer.c 
                       core/adios_bp_v1.c  
                       core/adios_endianness.c
                       core/futils.c 
                       core/adios_error.c 
                       core/bp_utils.c                
                       core/common_read.c 
                       core/globals.c 
                       core/adios_timing.c 
                       core/adios_read_hooks.c 
                       core/adios_transport_hooks.c 
                       core/adios_socket.c 
                       core/adios_logger.c 
                       core/util.c 
                       read/read_bp.c 
                       read/read_bp_staged.c 
                       read/read_bp_staged1.c 
                       write/adios_mpi.c 
                       write/adios_mpi_lustre.c 
                       write/adios_mpi_amr.c 
                       write/adios_posix.c 
                       write/adios_posix1.c)
  endif(BUILD_FORTRAN)

  set(libadios_a_CPPFLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")
  set(libadios_a_CFLAGS "${ADIOSLIB_CFLAGS}")
  set(libadios_nompi_a_CPPFLAGS "-D_NOMPI ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_SEQ_CPPFLAGS}") 
  set(libadios_nompi_a_CFLAGS "-D_NOMPI ${ADIOSLIB_SEQ_CFLAGS}")
  set(libadiosf_a_CPPFLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")
  set(libadiosf_a_CFLAGS "${ADIOSLIB_CFLAGS}")

  if(HAVE_DATASPACES)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_dataspaces.c read/read_dataspaces.c)
    set(libadios_nompi_a_SOURCES ${libadios_nompi_a_SOURCES} write/adios_dataspaces.c read/read_dataspaces.c)
    if(BUILD_FORTRAN)  
      set(FortranLibSources ${FortranLibSources} write/adios_dataspaces.c read/read_dataspaces.c)
    endif(BUILD_FORTRAN) 
  endif(HAVE_DATASPACES)

  if(HAVE_DIMES)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_dimes.c read/read_dimes.c)
    set(libadios_nompi_a_SOURCES ${libadios_nompi_a_SOURCES} write/adios_dimes.c read/read_dimes.c)
    if(BUILD_FORTRAN)  
      set(FortranLibSources ${FortranLibSources} write/adios_dimes.c read/read_dimes.c)
    endif(BUILD_FORTRAN)
  endif(HAVE_DIMES)

  if(HAVE_DATATAP)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_datatap.c read/read_datatap.c)
    set(libadios_nompi_a_SOURCES ${libadios_nompi_a_SOURCES} write/adios_datatap.c read/read_datatap.c)
    if(BUILD_FORTRAN)
      set(FortranLibSources ${FortranLibSources} write/adios_datatap.c read/read_datatap.c)
    endif(BUILD_FORTRAN)
  endif(HAVE_DATATAP)

  if(HAVE_PHDF5)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_phdf5.c)
    set(PHDF5_CPPFLAGS "${PHDF5_CPPFLAGS} -DPHDF5 -DH5_USE_16_API")
    set(libadios_a_CPPFLAGS "${libadios_a_CPPFLAGS} ${PHDF5_CPPFLAGS}")
    if(BUILD_FORTRAN)
      set(FortranLibSources ${FortranLibSources} write/adios_phdf5.c)
      set(PHDF5_CPPFLAGS "${PHDF5_CPPFLAGS} -DPHDF5 -DH5_USE_16_API")
      set(libadiosf_a_CPPFLAGS "${libadiosf_a_CPPFLAGS} ${PHDF5_CPPFLAGS}")
    endif(BUILD_FORTRAN) 
  endif(HAVE_PHDF5)

  if(HAVE_NC4PAR)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_nc4.c)
    if(BUILD_FORTRAN)
      set(FortranLibSources ${FortranLibSources} write/adios_nc4.c)
    endif(BUILD_FORTRAN)
  endif(HAVE_NC4PAR)

  if(HAVE_NSSI)
    set(libadios_a_SOURCES  ${libadios_a_SOURCES} nssi/adios_nssi_args.x)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} adios_nssi_args.c adios_nssi_args.h)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} nssi/adios_nssi_config.c write/adios_nssi.c nssi/adios_nssi_filter.c read/read_nssi.c)
    set(libadios_a_CPPFLAGS "${libadios_a_CPPFLAGS} -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL")
    if(BUILD_FORTRAN)
      set(libadiosf_a_SOURCES ${libadiosf_a_SOURCES} adios_nssi_args.c adios_nssi_args.h)
      set(libadiosf_a_SOURCES ${libadiosf_a_SOURCES} nssi/adios_nssi_args.x)
      set(FortranLibSources ${FortranLibSources} nssi/adios_nssi_config.c write/adios_nssi.c nssi/adios_nssi_filter.c read/read_nssi.c)
      set(libadiosf_a_CPPFLAGS ${libadiosf_a_CPPFLAGS} -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL)
    endif(BUILD_FORTRAN)

#not complete
#from line 200 to line 218
#not complete
    set(CLEANFILES ${CLEANFILES} adios_nssi_args.c adios_nssi_args.h)

  endif(HAVE_NSSI)


  if(RESEARCH_TRANSPORTS)
    set(libadios_a_SOURCES ${libadios_a_SOURCES}
                      write/adios_mpi_stagger.c 
                      write/adios_mpi_aggregate.c 
                      write/adios_provenance.c 
                      write/adios_mpi_cio.c 
                      write/adios_mpi_stripe.c 
                      write/adios_mpi_amr1.c 
                      write/adios_adaptive.ci)
    if(BUILD_FORTRAN)
      set(FortranLibSources ${FortranLibSources}
                          write/adios_mpi_stagger.c
                          write/adios_mpi_aggregate.c
                          write/adios_provenance.c
                          write/adios_mpi_cio.c
                          write/adios_mpi_stripe.c 
                          write/adios_mpi_amr1.c 
                          write/adios_adaptive.c)
    endif(BUILD_FORTRAN)
  else(RESEARCH_TRANSPORTS)
    set(libadios_a_CPPFLAGS "${libadios_a_CPPFLAGS} -DNO_RESEARCH_TRANSPORTS")
    set(libadios_nompi_a_CPPFLAGS "${libadios_nompi_a_CPPFLAGS} -DNO_RESEARCH_TRANSPORTS")
    set(libadiosf_a_CPPFLAGS "${libadiosf_a_CPPFLAGS} -DNO_RESEARCH_TRANSPORTS")
  endif(RESEARCH_TRANSPORTS)

  if(BUILD_FORTRAN)
    set(libadiosf_a_SOURCES ${FortranLibSources} core/adiosf_read.c)
    set(libadiosf_v1_a_SOURCES ${FortranLibSources} core/adiosf_read_v1.c)
    set(libadiosf_v1_a_CPPFLAGS "${libadiosf_a_CPPFLAGS}")
    set(libadiosf_v1_a_CFLAGS "${libadiosf_a_CFLAGS}")

    set(libadiosf_a_SOURCES ${libadiosf_a_SOURCES}
                         core/adiosf_defs_mod.F90
                         core/adiosf_write_mod.F90 
                         core/adiosf_read_mod.F90)

    set(nodist_include_HEADERS ${nodist_include_HEADERS} adios_write_mod.mod)
    set(CLEANFILES ${CLEANFILES} adios_write_mod.mod)
#not complete
#from line 261 to line 262
#not complet
    include_directories(${nodist_include_HEADERS})
    add_library(adiosf STATIC ${libadiosf_a_SOURCES})
    add_library(adiosf_v1 STATIC ${libadiosf_v1_a_SOURCES}) 
    set_target_properties(adiosf PROPERTIES COMPILE_FLAGS "${libadiosf_a_CPPFLAGS} ${libadiosf_a_CFLAGS}")
    set_target_properties(adiosf_v1 PROPERTIES COMPILE_FLAGS "${libadiosf_v1_a_CPPFLAGS} ${libadiosf_v1_a_CFLAGS}")
  endif(BUILD_FORTRAN)

  add_library(adios STATIC ${libadios_a_SOURCES})
  add_library(adios_nompi STATIC ${libadios_nompi_a_SOURCES})
  set_target_properties(adios PROPERTIES COMPILE_FLAGS "${libadiosf_a_CPPFLAGS} ${libadiosf_a_CFLAGS}") 
  set_target_properties(adios_nompi PROPERTIES COMPILE_FLAGS "${libadios_nompi_a_CPPFLAGS} ${libadios_nompi_a_CFLAGS}")
endif(BUILD_WRITE)
 

#start libadiosread.a
set(include_HEADERS ${include_HEADERS} public/adios.h 
                   public/adios_types.h 
                   public/adios_read.h 
                   public/adios_error.h 
                   public/mpidummy.h 
                   public/adios_read_v1.h 
                   public/adios_read_v1_defs.h 
                   public/adios_read_v2.h 
                   public/adios_selection.h) 

set(libadiosread_a_SOURCES core/adios_bp_v1.c
                      core/adios_endianness.c 
                      core/bp_utils.c 
                      core/futils.c 
                      core/adios_error.c 
                      core/adios_read.c 
                      core/adios_read_v1.c 
                      core/common_read.c 
                      core/globals.c 
                      core/adios_read_hooks.c 
                      core/adios_logger.c 
                      core/util.c 
                      read/read_bp.c 
                      read/read_bp_staged.c 
                      read/read_bp_staged1.c)

if(HAVE_DATASPACES)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} read/read_dataspaces.c)
endif(HAVE_DATASPACES)

if(HAVE_DIMES)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} read/read_dimes.c)
endif(HAVE_DIMES)

if(HAVE_DATATAP)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} read/read_datatap.c)
endif(HAVE_DATATAP)

if(HAVE_NSSI)
  set(dist_libadiosread_a_SOURCES nssi/adios_nssi_args.x)
  set(nodist_libadiosread_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} ${dist_libadiosread_a_SOURCES} ${nodist_libadiosread_a_SOURCES})
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} nssi/adios_nssi_config.c read/read_nssi.c)
  set(libadiosread_a_CPPFLAGS "${libadiosread_a_CPPFLAGS} -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL")
  set(libadiosread_a_CFLAGS "${libadiosread_a_CFLAGS} ${ADIOSLIB_CFLAGS}")
endif(HAVE_NSSI)

set(libadiosread_a_CPPFLAGS "${libadiosread_a_CPPFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS}")

#start libadiosreadf.a
if(BUILD_FORTRAN)
  set(FortranReadLibSource core/adios_bp_v1.c 
                      core/adios_endianness.c 
                      core/bp_utils.c 
                      core/futils.c 
                      core/adios_error.c 
                      core/common_read.c 
                      core/globals.c 
                      core/adios_read_hooks.c 
                      core/adios_logger.c 
                      core/util.c 
                      read/read_bp.c 
                      read/read_bp_staged.c 
                      read/read_bp_staged1.c)
  if(HAVE_DATASPACES)
    set(FortranReadLibSource ${FortranReadLibSource} read/read_dataspaces.c)
  endif(HAVE_DATASPACES)

  if(HAVE_DIMES)
    set(FortranReadLibSource ${FortranReadLibSource} read/read_dimes.c)
  endif(HAVE_DIMES)

  if(HAVE_DATATAP)
    set(FortranReadLibSource ${FortranReadLibSource} read/read_datatap.c)
  endif(HAVE_DATATAP)
  
  if(HAVE_NSSI)
    set(dist_libadiosreadf_a_SOURCES nssi/adios_nssi_args.x)
    set(nodist_libadiosreadf_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
    set(FortranReadLibSource ${FortranReadLibSource} ${dist_libadiosreadf_a_SOURCES} ${nodist_libadiosreadf_a_SOURCES})
    set(FortranReadLibSource ${FortranReadLibSource} nssi/adios_nssi_config.c read/read_nssi.c)
    set(libadiosreadf_a_CPPFLAGS "${libadiosreadf_a_CPPFLAGS} -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL")
    set(libadiosreadf_a_CFLAGS "${libadiosreadf_a_CFLAGS} ${ADIOSLIB_CFLAGS}")
  endif(HAVE_NSSI)

  set(libadiosreadf_a_CPPFLAGS "${libadiosreadf_a_CPPFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS}")
  set(libadiosreadf_a_SOURCES ${libadiosreadf_a_SOURCES} {FortranReadLibSource} core/adiosf_read.c)
  set(libadiosreadf_v1_a_CPPFLAGS "${libadiosreadf_a_CPPFLAGS}")
  set(libadiosreadf_v1_a_CFLAGS "${libadiosreadf_a_CFLAGS}")
  set(libadiosreadf_a_SOURCES ${libadiosreadf_a_SOURCES} core/adiosf_defs_mod.F90 core/adiosf_read_mod.F90)
  set(nodist_include_HEADERS ${nodist_include_HEADERS} adios_read_mod.mod adios_defs_mod.mod)
  set(CLEANFILES ${CLEANFILES} adios_read_mod.mod adios_defs_mod.mod)
  #not complete
  #from line 363 to line 367
  #not complete
  
endif(BUILD_FORTRAN)

#start libadiosread_nompi.a
set(libadiosread_nompi_a_SOURCES core/mpidummy.c
                      core/adios_bp_v1.c 
                      core/adios_endianness.c 
                      core/bp_utils.c 
                      core/futils.c 
                      core/adios_error.c 
                      core/adios_read.c 
                      core/adios_read_v1.c 
                      core/common_read.c 
                      core/adios_internals.c 
                      core/adios_logger.c 
                      core/buffer.c 
                      core/globals.c 
                      core/adios_read_hooks.c 
                      core/adios_transport_hooks.c 
                      core/util.c 
                      read/read_bp.c)

if(HAVE_DATASPACES)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} read/read_dataspaces.c)
endif(HAVE_DATASPACES)

if(HAVE_DIMES)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} read/read_dimes.c)
endif(HAVE_DIMES)

if(HAVE_DATATAP)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} read/read_datatap.c)
endif(HAVE_DATATAP)

if(HAVE_NSSI)
  set(dist_libadiosread_nompi_a_SOURCES nssi/adios_nssi_args.x)
  set(nodist_libadiosread_nompi_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} ${dist_libadiosread_nompi_a_SOURCES} ${nodist_libadiosread_nompi_a_SOURCES})
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} nssi/adios_nssi_config.c read/read_nssi.c)
  set(libadiosread_nompi_a_CPPFLAGS "${libadiosread_nompi_a_CPPFLAGS} -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NA")
endif(HAVE_NSSI)

set(libadiosread_nompi_a_CPPFLAGS "${libadiosread_nompi_a_CPPFLAGS} -D_NOMPI ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS}")

#start libadiosreadf_nompi.a
if(BUILD_FORTRAN)
  set(FortranReadSeqLibSource core/mpidummy.c
                          core/adios_bp_v1.c 
                          core/adios_endianness.c 
                          core/bp_utils.c 
                          core/futils.c 
                          core/adios_error.c 
                          core/adios_logger.c 
                          core/common_read.c 
                          core/globals.c 
                          core/adios_read_hooks.c 
                          core/util.c 
                          read/read_bp.c)
  if(HAVE_DATASPACES)
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} read/read_dataspaces.c)
  endif(HAVE_DATASPACES)

  if(HAVE_DIMES)
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} read/read_dimes.c)
  endif(HAVE_DIMES)

  if(HAVE_DATATAP)
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} read/read_datatap.c)
  endif(HAVE_DATATAP)

  if(HAVE_NSSI)
    set(dist_libadiosreadf_nompi_a_SOURCES nssi/adios_nssi_args.x)
    set(nodist_libadiosreadf_nompi_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} ${dist_libadiosreadf_nompi_a_SOURCES} ${nodist_libadiosreadf_nompi_a_SOURCES})
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} nssi/adios_nssi_config.c read/read_nssi.c)
    set(libadiosreadf_nompi_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS} -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL")
    set(libadiosreadf_nompi_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")
    set(libadiosreadf_nompi_a_CFLAGS "${ADIOSLIB_CFLAGS}")
  endif(HAVE_NSSI)

  set(libadiosreadf_nompi_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS} -D_NOMPI ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS}")
  set(libadiosreadf_nompi_a_SOURCES ${FortranReadSeqLibSource} core/adiosf_read.c)
  set(libadiosreadf_nompi_v1_a_SOURCES ${FortranReadSeqLibSource} core/adiosf_read_v1.c)
  set(libadiosreadf_nompi_v1_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS}")
  set(libadiosreadf_nompi_v1_a_CFLAGS "${libadiosreadf_nompi_a_CFLAGS}")

endif(BUILD_FORTRAN)

#start libadios_internal_nompi.a
set(libadios_internal_nompi_a_SOURCES core/mpidummy.c 
                                    core/adios_bp_v1.c 
                                    core/adios_endianness.c 
                                    core/bp_utils.c 
                                    core/adios_internals.c 
                                    core/buffer.c 
                                    core/adios_error.c 
                                    core/adios_logger.c 
                                    core/adios_internals_mxml.c 
                                    core/adios_timing.c 
                                    core/util.c 
                                    core/futils.c 
                                    core/adios_transport_hooks.c)

set(libadios_internal_nompi_a_CPPFLAGS "-D_INTERNAL ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")

if(HAVE_DMALLOC)
  set(libadiosread_a_CPPFLAGS "${libadiosread_a_CPPFLAGS} -DDMALLOC")
  set(libadiosread_nompi_a_CPPFLAGS "${libadiosread_nompi_a_CPPFLAGS} -DDMALLOC") 
endif(HAVE_DMALLOC)

#include_directories(${include_HEADERS})
add_library(adiosread STATIC ${libadiosread_a_SOURCES})
set_target_properties(adiosread PROPERTIES COMPILE_FLAGS "${libadiosread_a_CPPFLAGS} ${libadiosread_a_CFLAGS}")

add_library(adiosread_nompi ${libadiosread_nompi_a_SOURCES})
set_target_properties(adiosread_nompi PROPERTIES COMPILE_FLAGS "${libadiosread_nompi_a_CPPFLAGS} ${libadiosread_nompi_a_CFLAGS}")

add_library(adios_internal_nompi ${libadios_internal_nompi_a_SOURCES})
set_target_properties(adios_internal_nompi PROPERTIES COMPILE_FLAGS "${libadios_internal_nompi_a_CPPFLAGS}")

if(BUILD_FORTRAN)
  include_directories(${nodist_include_HEADERS})
  add_library(adiosreadf ${libadiosreadf_a_SOURCES})
  set_target_properties(adiosreadf PROPERTIES COMPILE_FLAGS "${libadiosreadf_a_CPPFLAGS} ${libadiosreadf_a_CFLAGS}")

  add_library(adiosreadf_v1 ${libadiosreadf_v1_a_SOURCES} ${nodist_include_HEADERS})
  set_target_properties(libadiosreadf_v1 PROPERTIES COMPILE_FLAGS "${libadiosreadf_v1_a_CPPFLAGS} ${libadiosreadf_v1_a_CFLAGS}")

  add_library(adiosreadf_nompi ${libadiosreadf_nompi_a_SOURCES})
  set_target_properties(adiosreadf_nompi  PROPERTIES COMPILE_FLAGS "${libadiosreadf_nompi_a_CPPFLAGS} ${libadiosreadf_nompi_a_CFLAGS}")

  add_library(adiosreadf_nompi_v1 ${libadiosreadf_nompi_v1_a_SOURCES})
  set_target_properties(adiosreadf_nompi_v1 PROPERTIES COMPILE_FLAGS "${libadiosreadf_nompi_v1_a_CPPFLAGS} ${libadiosreadf_nompi_v1_a_CFLAGS}")

endif(BUILD_FORTRAN)

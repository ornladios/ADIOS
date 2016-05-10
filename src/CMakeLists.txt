include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/src/core ${PROJECT_SOURCE_DIR}/src/core/transforms)
include_directories(${PROJECT_BINARY_DIR}/src ${PROJECT_BINARY_DIR}/src/public)
link_directories(${PROJECT_BINARY_DIR}/src)

if(HAVE_NSSI)
  add_executable(nssi-staging-server nssi/nssi_staging_server.cpp nssi/aggregation.cpp)
  set_target_properties(nssi-staging-server PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NAL") 
  target_link_libraries(nssi-staging-server adios.a ${ADIOSLIB_LDADD} ${NSSI_SERVER_LIBS})

  add_executable(nssi-coupling-server nssi/nssi_coupling_server.cpp)
  set_target_properties(nssi-coupling-server PROPERTIES COMPILE_FLAGS "${ADIOSLIB_CFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NAL") 
  target_link_libraries(nssi-coupling-server adios.a ${ADIOSLIB_LDADD} ${NSSI_SERVER_LIBS})
endif(HAVE_NSSI)



set(CLEANFILES "")
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

## Transform plugin source files
set(transforms_write_method_SOURCES "")
set(transforms_read_method_SOURCES "")

#include (transforms/Makefile.plugins.cmake)

##
## Core transform framework source files
##
set (transforms_common_HDRS core/adios_copyspec.h 
  core/adios_subvolume.h 
  core/adios_selection_util.h 
  core/transforms/adios_transforms_common.h 
  core/transforms/adios_transforms_hooks.h 
  core/transforms/adios_transforms_util.h 
  core/adios_subvolume.h)

set (transforms_read_HDRS core/transforms/adios_transforms_read.h 
  core/transforms/adios_transforms_hooks_read.h 
  core/transforms/adios_transforms_reqgroup.h 
  core/transforms/adios_transforms_datablock.h
  core/transforms/adios_transforms_transinfo.h
  core/transforms/adios_patchdata.h)

set (transforms_write_HDRS  core/transforms/adios_transforms_write.h 
  core/transforms/adios_transforms_hooks_write.h 
  core/transforms/adios_transforms_specparse.h)

set (transforms_common_SOURCES  ${transforms_common_HDRS} 
  core/transforms/adios_transforms_common.c 
  core/transforms/adios_transforms_hooks.c 
  core/adios_copyspec.c 
  core/adios_subvolume.c 
  core/transforms/plugindetect/detect_plugin_infos.h 
  core/transforms/plugindetect/detect_plugin_types.h 
  core/transforms/plugindetect/plugin_info_types.h)

set (transforms_read_SOURCES ${transforms_read_HDRS} 
  core/transforms/adios_transforms_read.c 
  core/transforms/adios_transforms_hooks_read.c 
  core/transforms/adios_transforms_reqgroup.c 
  core/transforms/adios_transforms_datablock.c 
  core/transforms/adios_patchdata.c 
  transforms/adios_transform_alacrity_read.c
  transforms/adios_transform_isobar_read.c
  transforms/adios_transform_aplod_read.c
  transforms/adios_transform_szip_read.c 
  transforms/adios_transform_bzip2_read.c
  transforms/adios_transform_identity_read.c
  transforms/adios_transform_zlib_read.c
  core/adios_selection_util.c 
  core/transforms/plugindetect/detect_plugin_read_hook_decls.h
  core/transforms/plugindetect/detect_plugin_read_hook_reg.h
  core/adios_read_ext.c
  ${transforms_read_method_SOURCES})

set (transforms_write_SOURCES ${transforms_write_HDRS} 
  core/transforms/adios_transforms_write.c 
  core/transforms/adios_transforms_hooks_write.c 
  core/transforms/adios_transforms_util.c 
  core/transforms/adios_transforms_specparse.c 
  transforms/adios_transform_alacrity_write.c
  transforms/adios_transform_aplod_write.c
  transforms/adios_transform_bzip2_write.c
  transforms/adios_transform_identity_write.c
  transforms/adios_transform_isobar_write.c
  transforms/adios_transform_szip_write.c
  transforms/adios_transform_zlib_write.c
  ${transforms_write_method_SOURCES})

#######Query source files
set(query_common_HDRS query/common_query.h
  query/adios_query_hooks.h
  query/query_utils.h)

set(query_common_SOURCES ${query_common_HDRS}
  query/common_query.c
  query/adios_query_hooks.c
  query/query_utils.c)

# Include source files that are specific to each query plugin
set(query_method_HDRS "")
set(query_method_SOURCES "")
include(query/Makefile.plugins.cmake)

set(query_method_SOURCES ${query_method_SOURCES} ${query_method_HDRS})

set(query_SOURCES ${query_common_SOURCES} ${query_method_SOURCES})

set(query_C_SOURCES ${query_SOURCES} query/adios_query.c)
set(query_F_SOURCES ${query_SOURCES} query/adiosf_query.c)

#start adios.a and adios_nompi.a
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
    core/adios_infocache.c
    core/adios_read_ext.c
    core/globals.c 
    core/adios_timing.c 
    core/adios_read_hooks.c 
    core/adios_transport_hooks.c 
    core/adios_socket.c 
    core/adios_logger.c 
    core/qhashtbl.c 
    ${transforms_common_SOURCES} 
    ${transforms_read_SOURCES} 
    ${transforms_write_SOURCES} 
    ${query_C_SOURCES}
    core/util.c core/strutil.c 
    core/qhashtbl.c 
    read/read_bp.c 
    read/read_bp_staged.c 
    read/read_bp_staged1.c
    write/adios_mpi.c
    write/adios_mpi_lustre.c
    write/adios_mpi_amr.c
    write/adios_posix.c
    write/adios_var_merge.c write/adios_sirius.c)

  if(HAVE_BGQ)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_mpi_bgq.c)
  endif()

  set(libadios_nompi_a_SOURCES core/adios.c 
    core/common_adios.c 
    core/adios_internals.c 
    core/adios_internals_mxml.c 
    ${transforms_common_SOURCES} 
    ${transforms_read_SOURCES} 
    ${transforms_write_SOURCES} 
    ${query_C_SOURCES}
    core/buffer.c 
    core/adios_bp_v1.c  
    core/adios_endianness.c 
    core/bp_utils.c 
    core/futils.c 
    core/adios_error.c 
    core/adios_read.c 
    core/adios_read_v1.c 
    core/common_read.c 
    core/adios_infocache.c
    core/adios_read_ext.c
    core/globals.c 
    core/mpidummy.c 
    core/adios_timing.c 
    core/adios_read_hooks.c 
    core/adios_transport_hooks.c 
    core/adios_socket.c 
    core/adios_logger.c 
    core/util.c core/strutil.c 
    core/qhashtbl.c 
    read/read_bp.c 
    read/read_bp_staged.c 
    read/read_bp_staged1.c 
    write/adios_posix.c)

  #start adiosf.a and adiosf_v1.a
  if(BUILD_FORTRAN)
    set(FortranLibSources core/adiosf.c 
      core/common_adios.c 
      core/adios_internals.c 
      core/adios_internals_mxml.c
      ${transforms_common_SOURCES} 
      ${transforms_read_SOURCES} 
      ${transforms_write_SOURCES} 
      ${query_F_SOURCES}
      core/buffer.c 
      core/adios_bp_v1.c  
      core/adios_endianness.c
      core/futils.c 
      core/adios_error.c 
      core/bp_utils.c                
      core/common_read.c 
      core/adios_infocache.c
      core/adios_read_ext.c
      core/globals.c 
      core/adios_timing.c 
      core/adios_read_hooks.c 
      core/adios_transport_hooks.c 
      core/adios_socket.c 
      core/adios_logger.c 
      core/util.c core/strutil.c 
      core/qhashtbl.c 
      read/read_bp.c 
      read/read_bp_staged.c 
      read/read_bp_staged1.c 
      write/adios_posix.c)

    set(FortranLibMPISources write/adios_mpi.c
      write/adios_mpi_lustre.c
      write/adios_mpi_amr.c
      write/adios_var_merge.c write/adios_sirius.c)
    if(HAVE_BGQ)
      set(FortranLibMPISources ${FortranLibMPISources} write/adios_mpi_bgq.c)
    endif(HAVE_BGQ)
    set(FortranLibSEQSources core/mpidummy.c)
  endif(BUILD_FORTRAN)

  set(libadios_a_CPPFLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")
  set(libadios_a_CFLAGS "${ADIOSLIB_CFLAGS}")
  set(libadios_nompi_a_CPPFLAGS "${MACRODEFFLAG}_NOMPI ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_SEQ_CPPFLAGS}") 
  set(libadios_nompi_a_CFLAGS "${MACRODEFFLAG}_NOMPI ${ADIOSLIB_SEQ_CFLAGS}")
  set(libadiosf_a_CPPFLAGS "${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")
  set(libadiosf_a_CFLAGS "${ADIOSLIB_CFLAGS}")
  set(libadiosf_nompi_a_CPPFLAGS "${MACRODEFFLAG}_NOMPI ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_SEQ_CPPFLAGS}")
  set(libadiosf_nompi_a_CFLAGS "${MACRODEFFLAG}_NOMPI ${ADIOSLIB_SEQ_CFLAGS}")

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
    set(PHDF5_FLAGS "${PHDF5_FLAGS} ${MACRODEFFLAG}PHDF5 ${MACRODEFFLAG}H5_USE_16_API")
    set(libadios_a_CPPFLAGS "${libadios_a_CPPFLAGS} ${PHDF5_FLAGS}")
    if(BUILD_FORTRAN)
      set(FortranLibMPISources ${FortranLibMPISources} write/adios_phdf5.c)
      set(libadiosf_a_CPPFLAGS "${libadiosf_a_CPPFLAGS} ${PHDF5_FLAGS}")
    endif(BUILD_FORTRAN) 
  endif(HAVE_PHDF5)

  if(HAVE_NC4PAR)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_nc4.c)
    if(BUILD_FORTRAN)
      set(FortranLibMPISources ${FortranLibMPISources} write/adios_nc4.c)
    endif(BUILD_FORTRAN)
  endif(HAVE_NC4PAR)

  if(HAVE_FLEXPATH)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_flexpath.c read/read_flexpath.c)
    set(libadios_nompi_a_SOURCES ${libadios_nompi_a_SOURCES} write/adios_flexpath.c read/read_flexpath.c)
    if(BUILD_FORTRAN)
      set(FortranLibSources ${FortranLibSources} write/adios_flexpath.c read/read_flexpath.c)
    endif(BUILD_FORTRAN)
  endif(HAVE_FLEXPATH)

  if(HAVE_ICEE)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} write/adios_icee.c read/read_icee.c)
    set(libadios_nompi_a_SOURCES ${libadios_nompi_a_SOURCES} write/adios_icee.c  read/read_icee.c)
    if(BUILD_FORTRAN)
      set(FortranLibSources ${FortranLibSources} write/adios_icee.c read/read_icee.c)
    endif()
  endif()

  if(HAVE_NSSI)
    set(dist_libadios_a_SOURCES  nssi/adios_nssi_args.x)
    set(nodist_libadios_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
    set(libadios_a_SOURCES ${libadios_a_SOURCES} ${dist_libadios_a_SOURCES} ${nodist_libadios_a_SOURCES})
    set(libadios_a_SOURCES ${libadios_a_SOURCES} nssi/adios_nssi_config.c write/adios_nssi.c nssi/adios_nssi_filter.c read/read_nssi.c)
    set(libadios_a_CPPFLAGS "${libadios_a_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NAL")
    if(BUILD_FORTRAN)
      set(nodist_libadiosf_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
      set(dist_libadiosf_a_SOURCES nssi/adios_nssi_args.x)
      set(libadiosf_a_SOURCES ${libadiosf_a_SOURCES} ${nodist_libadiosf_a_SOURCES} ${dist_libadiosf_a_SOURCES})
      set(FortranLibSources ${FortranLibSources} nssi/adios_nssi_config.c write/adios_nssi.c nssi/adios_nssi_filter.c read/read_nssi.c)
      set(libadiosf_a_CPPFLAGS "${libadiosf_a_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NAL")
    endif(BUILD_FORTRAN)

    #adios_nssi_config.c: adios_nssi_args.c
    #use include_directories(${PROJECT_BINARY_DIR}/src) to let adios_nssi_config.c
    #find the header needed

    add_custom_command (
      OUTPUT adios_nssi_args.h
      COMMAND rm -f $@
      COMMAND rpcgen -h ${PROJECT_SOURCE_DIR}/src/nssi/adios_nssi_args.x >> $@
      DEPENDS ${PROJECT_SOURCE_DIR}/src/nssi/adios_nssi_args.x ${PROJECT_SOURCE_DIR}/src/CMakeLists.txt
      )

    add_custom_command (
      OUTPUT adios_nssi_args.c
      COMMAND rm -f adios_nssi_args.c
      COMMAND rpcgen -c ${PROJECT_SOURCE_DIR}/src/nssi/adios_nssi_args.x | grep -v "include \"\.\." > adios_nssi_args.c
      DEPENDS ${PROJECT_SOURCE_DIR}/src/nssi/adios_nssi_args.h
      )
    
    set(CLEANFILES ${CLEANFILES} adios_nssi_args.c adios_nssi_args.h)

  endif(HAVE_NSSI)


  if(RESEARCH_TRANSPORTS)
    #        set(libadios_a_SOURCES ${libadios_a_SOURCES}
    #                      write/adios_mpi_stagger.c 
    #                      write/adios_mpi_aggregate.c 
    #                      write/adios_provenance.c 
    #                      write/adios_mpi_cio.c 
    #                      write/adios_mpi_stripe.c 
    #                      write/adios_mpi_amr1.c 
    #                      write/adios_adaptive.c)
    if(BUILD_FORTRAN)
      #            set(FortranLibMPISources ${FortranLibMPISources}
      #                          write/adios_mpi_stagger.c
      #                          write/adios_mpi_aggregate.c
      #                          write/adios_provenance.c
      #                          write/adios_mpi_cio.c
      #                          write/adios_mpi_stripe.c 
      #                          write/adios_mpi_amr1.c 
      #                          write/adios_adaptive.c)
    endif(BUILD_FORTRAN)
  else(RESEARCH_TRANSPORTS)
    set(libadios_a_CPPFLAGS "${libadios_a_CPPFLAGS} ${MACRODEFFLAG}NO_RESEARCH_TRANSPORTS")
    set(libadios_nompi_a_CPPFLAGS "${libadios_nompi_a_CPPFLAGS} ${MACRODEFFLAG}NO_RESEARCH_TRANSPORTS")
    set(libadiosf_a_CPPFLAGS "${libadiosf_a_CPPFLAGS} ${MACRODEFFLAG}NO_RESEARCH_TRANSPORTS")
    set(libadiosf_nompi_a_CPPFLAGS "${libadiosf_nompi_a_CPPFLAGS} ${MACRODEFFLAG}NO_RESEARCH_TRANSPORTS")
  endif(RESEARCH_TRANSPORTS)

  # add library adios and adios_nompi
  add_library(adios STATIC ${libadios_a_SOURCES})
  add_library(adios_nompi STATIC ${libadios_nompi_a_SOURCES})
  set_target_properties(adios PROPERTIES COMPILE_FLAGS "${libadios_a_CPPFLAGS} ${libadios_a_CFLAGS}")
  set_target_properties(adios_nompi PROPERTIES COMPILE_FLAGS "${libadios_nompi_a_CPPFLAGS} ${libadios_nompi_a_CFLAGS}")
  #   message(STATUS "${libadios_a_CPPFLAGS}")
  #   message(STATUS "${libadios_a_CFLAGS}")

  # Build four different Fortran libraries, for the two read APIs
  # and their sequential versions
  if(BUILD_FORTRAN)
    set(libadiosf_a_SOURCES ${libadiosf_a_SOURCES} ${FortranLibSources} ${FortranLibMPISources} core/adiosf_read.c)
    set(libadiosf_v1_a_SOURCES ${FortranLibSources} ${FortranLibMPISources} core/adiosf_read_v1.c)
    set(libadiosf_nompi_a_SOURCES ${FortranLibSources} ${FortranLibSEQSources} core/adiosf_read.c)
    set(libadiosf_nompi_v1_a_SOURCES ${FortranLibSources} ${FortranLibSEQSources} core/adiosf_read_v1.c)
    set(libadiosf_v1_a_CPPFLAGS "${libadiosf_a_CPPFLAGS}")
    set(libadiosf_v1_a_CFLAGS "${libadiosf_a_CFLAGS}")
    set(libadiosf_nompi_v1_a_CPPFLAGS "${libadiosf_nompi_a_CPPFLAGS}")
    set(libadiosf_nompi_v1_a_CFLAGS "${libadiosf_nompi_a_CFLAGS}")

    set(libadiosf_a_SOURCES ${libadiosf_a_SOURCES}
      core/adiosf_defs_mod.f90
      core/adiosf_write_mod.f90 
      core/adiosf_read_mod.f90
      query/adios_query_mod.f90)
    
    set(libadiosf_nompi_a_SOURCES ${libadiosf_nompi_a_SOURCES}
      core/adiosf_defs_mod.f90
      core/adiosf_write_mod.f90
      core/adiosf_read_mod.f90
      query/adios_query_mod.f90)

    set(libadiosf_v1_a_SOURCES ${libadiosf_v1_a_SOURCES} 
      core/adiosf_defs_mod.f90
      core/adiosf_write_mod.f90)

    set(libadiosf_nompi_v1_a_SOURCES ${libadiosf_nompi_v1_a_SOURCES}
      core/adiosf_defs_mod.f90
      core/adiosf_write_mod.f90)

    set(nodist_include_HEADERS ${nodist_include_HEADERS} adios_write_mod.mod)
    set(CLEANFILES ${CLEANFILES} adios_write_mod.mod)

    #       add_custom_command(
    #           OUTPUT adios_write_mod.mod
    #           COMMAND $(FC) -c core/adiosf_write_mod.f90
    #           DEPENDS core/adiosf_write_mod.f90 adios_defs_mod.mod
    #           )
    #       add_custom_command(
    #           OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/adios_write_mod.mod
    #           COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/adios_write_mod.mod ${CMAKE_CURRENT_SOURCE_DIR}
    #           DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/adios_write_mod.mod
    #           )

    # add_custom_command
    # (
    #             OUTPUT adios_write_mod.o
    # #            COMMAND  $(FC) -c core/adiosf_write_mod.f90
    #             DEPENDS adios_defs_mod.o
    #         )

    
    # add_custom_command(
    #     OUTPUT adios_write_mod.mod
    #     DEPENDS adios_defs_mod.o
    # )

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/adios_write_mod.mod
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/adios_write_mod.mod ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS adios_write_mod.mod
      )
    #       file(COPY ${CMAKE_CURRENT_BINARY_DIR}/adios_write_mod.mod DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

    ##include_directories(${nodist_include_HEADERS})
    #       install(FILES ${nodist_include_HEADERS} DESTINATION ${includedir})
    install(FILES ${PROJECT_BINARY_DIR}/src/adios_write_mod.mod ${PROJECT_BINARY_DIR}/src/adios_defs_mod.mod DESTINATION ${includedir})
    ##add library adiosf and adiosf_v1
    add_library(adiosf STATIC ${libadiosf_a_SOURCES})
    add_library(adiosf_v1 STATIC ${libadiosf_v1_a_SOURCES}) 
    add_library(adiosf_nompi STATIC ${libadiosf_nompi_a_SOURCES})
    add_library(adiosf_nompi_v1 STATIC ${libadiosf_nompi_v1_a_SOURCES})
    set_target_properties(adiosf PROPERTIES COMPILE_FLAGS "${libadiosf_a_CPPFLAGS} ${libadiosf_a_CFLAGS}")
    set_target_properties(adiosf_v1 PROPERTIES COMPILE_FLAGS "${libadiosf_v1_a_CPPFLAGS} ${libadiosf_v1_a_CFLAGS}")
    set_target_properties(adiosf_nompi PROPERTIES COMPILE_FLAGS "${libadiosf_nompi_a_CPPFLAGS} ${libadiosf_nompi_a_CFLAGS}")
    set_target_properties(adiosf_nompi_v1 PROPERTIES COMPILE_FLAGS "${libadiosf_nompi_v1_a_CPPFLAGS} ${libadiosf_nompi_v1_a_CFLAGS}")
    #       message(STATUS "${libadiosf_a_CPPFLAGS}")
    #       message(STATUS "${libadiosf_a_CFLAGS}")
  endif(BUILD_FORTRAN)

endif(BUILD_WRITE)


#start libadiosread.a
set(include_HEADERS ${include_HEADERS} public/adios.h 
  public/adios_types.h 
  public/adios_read.h 
  public/adios_error.h
  public/adios_mpi.h 
  public/mpidummy.h 
  public/adios_read_v1.h 
  public/adios_read_v1_defs.h 
  public/adios_read_v2.h 
  public/adios_read_v2_fwd.h 
  public/adios_selection.h
  public/adios_schema.h
  public/adios_link.h
  public/adios_read_ext.h
  public/adios_query.h)

set(libadiosread_a_SOURCES core/adios_bp_v1.c
  core/adios_endianness.c 
  core/bp_utils.c 
  core/futils.c 
  core/adios_error.c 
  core/adios_read.c 
  core/adios_read_v1.c 
  core/common_read.c 
  core/adios_infocache.c
  core/adios_read_ext.c
  ${transforms_common_SOURCES} 
  ${transforms_read_SOURCES} 
  ${query_C_SOURCES}
  core/globals.c 
  core/adios_read_hooks.c 
  core/adios_logger.c 
  core/util.c core/strutil.c 
  core/qhashtbl.c 
  read/read_bp.c 
  read/read_bp_staged.c 
  read/read_bp_staged1.c)

if(HAVE_DMALLOC)
  set(libadiosread_a_CPPFLAGS "${libadiosread_a_CPPFLAGS} ${MACRODEFFLAG}DMALLOC")
endif(HAVE_DMALLOC)

if(HAVE_DATASPACES)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} read/read_dataspaces.c)
endif(HAVE_DATASPACES)

if(HAVE_DIMES)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} read/read_dimes.c)
endif(HAVE_DIMES)

if(HAVE_DATATAP)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} read/read_datatap.c)
endif(HAVE_DATATAP)

if(HAVE_FLEXPATH)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} read/read_flexpath.c)
endif(HAVE_FLEXPATH)

if(HAVE_ICEE)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} read/read_icee.c)
endif()

if(HAVE_NSSI)
  set(dist_libadiosread_a_SOURCES nssi/adios_nssi_args.x)
  set(nodist_libadiosread_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} ${dist_libadiosread_a_SOURCES} ${nodist_libadiosread_a_SOURCES})
  set(libadiosread_a_SOURCES ${libadiosread_a_SOURCES} nssi/adios_nssi_config.c read/read_nssi.c)
  set(libadiosread_a_CPPFLAGS "${libadiosread_a_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NAL")
  set(libadiosread_a_CFLAGS "${libadiosread_a_CFLAGS} ${ADIOSLIB_CFLAGS}")
endif(HAVE_NSSI)

set(libadiosread_a_CPPFLAGS "${libadiosread_a_CPPFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS}")

#add library adiosread
#include_directories(${include_HEADERS})
install(FILES ${include_HEADERS} DESTINATION ${includedir})
# install generated header file adios_version.h
install(FILES ${PROJECT_BINARY_DIR}/src/public/adios_version.h DESTINATION ${includedir})
add_library(adiosread STATIC ${libadiosread_a_SOURCES})
set_target_properties(adiosread PROPERTIES COMPILE_FLAGS "${libadiosread_a_CPPFLAGS} ${libadiosread_a_CFLAGS}")

#start libadiosreadf.a libadiosreadf_v1.a
if(BUILD_FORTRAN)
  set(FortranReadLibSource core/adios_bp_v1.c 
    core/adios_endianness.c 
    core/bp_utils.c 
    core/futils.c 
    core/adios_error.c 
    core/common_read.c 
    core/adios_infocache.c
    core/adios_read_ext.c
    ${transforms_common_SOURCES} 
    ${transforms_read_SOURCES} 
    ${query_F_SOURCES}
    core/globals.c 
    core/adios_read_hooks.c 
    core/adios_logger.c 
    core/util.c core/strutil.c 
    core/qhashtbl.c 
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

  if(HAVE_FLEXPATH)
    set(FortranReadLibSource ${FortranReadLibSource} read/read_flexpath.c)
  endif(HAVE_FLEXPATH)
  
  if(HAVE_ICEE)
    set(FortranReadLibSource ${FortranReadLibSource} read/read_icee.c)
  endif()

  if(HAVE_NSSI)
    set(dist_libadiosreadf_a_SOURCES nssi/adios_nssi_args.x)
    set(nodist_libadiosreadf_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
    set(libadiosreadf_a_SOURCES ${libadiosreadf_a_SOURCES} ${dist_libadiosreadf_a_SOURCES} ${nodist_libadiosreadf_a_SOURCES})
    set(FortranReadLibSource ${FortranReadLibSource} nssi/adios_nssi_config.c read/read_nssi.c)
    set(libadiosreadf_a_CPPFLAGS "${libadiosreadf_a_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NAL")
    set(libadiosreadf_a_CFLAGS "${libadiosreadf_a_CFLAGS} ${ADIOSLIB_CFLAGS}")
  endif(HAVE_NSSI)

  set(libadiosreadf_a_CPPFLAGS "${libadiosreadf_a_CPPFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS}")
  set(libadiosreadf_a_SOURCES ${libadiosreadf_a_SOURCES} ${FortranReadLibSource} core/adiosf_read.c)
  set(libadiosreadf_v1_a_SOURCES ${FortranReadLibSource} core/adiosf_read_v1.c)
  set(libadiosreadf_v1_a_CPPFLAGS "${libadiosreadf_a_CPPFLAGS}")
  set(libadiosreadf_v1_a_CFLAGS "${libadiosreadf_a_CFLAGS}")
  set(libadiosreadf_a_SOURCES ${libadiosreadf_a_SOURCES}
    core/adiosf_defs_mod.f90 
    core/adiosf_read_mod.f90
    query/adios_query_mod.f90)
  set(nodist_include_HEADERS ${nodist_include_HEADERS} adios_read_mod.mod adios_defs_mod.mod adios_query_mod.mod)
  set(CLEANFILES ${CLEANFILES} adios_read_mod.mod adios_defs_mod.mod adios_query_mod.mod)
  
  # add_custom_command(
  #         OUTPUT adios_read_mod.o
  # #        COMMAND $(FC) -c core/adiosf_read_mod.f90
  #         DEPENDS adios_defs_mod.o
  #         )  
  #     add_custom_command(
  #         OUTPUT adios_read_mod.mod
  #         DEPENDS adios_read_mod.o
  # )
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/adios_read_mod.mod
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/adios_read_mod.mod ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS adios_read_mod.mod
    )

  # add_custom_command(
  #     OUTPUT adios_query_mod.o
  #     DEPENDS adios_read_mod.o
  #     )
  # add_custom_command(
  #     OUTPUT adios_query_mod.mod
  #     DEPENDS adios_query_mod.o
  #)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/adios_query_mod.mod
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/adios_query_mod.mod ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS adios_query_mod.mod
    ) 

  # add_custom_command(
  #     OUTPUT adios_defs_mod.mod
  #     DEPENDS adios_defs_mod.o
  # )
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/adios_defs_mod.mod
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/adios_defs_mod.mod ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS adios_defs_mod.mod
    )
  #   file(COPY ${CMAKE_CURRENT_BINARY_DIR}/adios_defs_mod.mod DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
  #add library adiosreadf and adiosreadf_v1
  #include_directories(${nodist_include_HEADERS})
  #  install(FILES ${nodist_include_HEADERS} DESTINATION ${includedir})
  install(FILES ${PROJECT_BINARY_DIR}/src/adios_read_mod.mod DESTINATION ${includedir})
  add_library(adiosreadf ${libadiosreadf_a_SOURCES})
  set_target_properties(adiosreadf PROPERTIES COMPILE_FLAGS "${libadiosreadf_a_CPPFLAGS} ${libadiosreadf_a_CFLAGS}")
  
  #  message(STATUS "${nodist_include_HEADERS}")
  #  message(STATUS "${libadiosreadf_v1_a_SOURCES}")
  add_library(adiosreadf_v1 ${libadiosreadf_v1_a_SOURCES})
  set_target_properties(adiosreadf_v1 PROPERTIES COMPILE_FLAGS "${libadiosreadf_v1_a_CPPFLAGS} ${libadiosreadf_v1_a_CFLAGS}")

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
  core/adios_infocache.c
  core/adios_read_ext.c
  ${transforms_common_SOURCES} 
  ${transforms_read_SOURCES} 
  ${query_C_SOURCES}
  #                      core/adios_internals.c 
  core/adios_logger.c 
  core/buffer.c 
  core/globals.c 
  core/adios_read_hooks.c 
  #                      core/adios_transport_hooks.c 
  core/util.c core/strutil.c 
  core/qhashtbl.c 
  read/read_bp.c)

if(HAVE_DMALLOC)
  set(libadiosread_nompi_a_CPPFLAGS "${libadiosread_nompi_a_CPPFLAGS} ${MACRODEFFLAG}DMALLOC")
endif(HAVE_DMALLOC)

if(HAVE_DATASPACES)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} read/read_dataspaces.c)
endif(HAVE_DATASPACES)

if(HAVE_DIMES)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} read/read_dimes.c)
endif(HAVE_DIMES)

if(HAVE_DATATAP)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} read/read_datatap.c)
endif(HAVE_DATATAP)

if(HAVE_FLEXPATH)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} read/read_flexpath.c)
endif(HAVE_FLEXPATH)

if(HAVE_ICEE)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} read/read_icee.c)
endif()

if(HAVE_NSSI)
  set(dist_libadiosread_nompi_a_SOURCES nssi/adios_nssi_args.x)
  set(nodist_libadiosread_nompi_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} ${dist_libadiosread_nompi_a_SOURCES} ${nodist_libadiosread_nompi_a_SOURCES})
  set(libadiosread_nompi_a_SOURCES ${libadiosread_nompi_a_SOURCES} nssi/adios_nssi_config.c read/read_nssi.c)
  set(libadiosread_nompi_a_CPPFLAGS "${libadiosread_nompi_a_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NA")
endif(HAVE_NSSI)

set(libadiosread_nompi_a_CPPFLAGS "${libadiosread_nompi_a_CPPFLAGS} ${MACRODEFFLAG}_NOMPI ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS}")

#add library adiosread_nompi
add_library(adiosread_nompi ${libadiosread_nompi_a_SOURCES})
set_target_properties(adiosread_nompi PROPERTIES COMPILE_FLAGS "${libadiosread_nompi_a_CPPFLAGS} ${libadiosread_nompi_a_CFLAGS}")


#start libadiosreadf_nompi.a libadiosreadf_nompi_v1.a
if(BUILD_FORTRAN)
  set(FortranReadSeqLibSource core/mpidummy.c
    core/adios_bp_v1.c 
    core/adios_endianness.c 
    core/bp_utils.c 
    core/futils.c 
    core/adios_error.c 
    core/adios_logger.c 
    core/common_read.c 
    core/adios_infocache.c
    core/adios_read_ext.c
    ${transforms_common_SOURCES} 
    ${transforms_read_SOURCES} 
    ${query_F_SOURCES}
    core/globals.c 
    core/adios_read_hooks.c 
    core/util.c core/strutil.c 
    core/qhashtbl.c 
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

  if(HAVE_FLEXPATH)
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} read/read_flexpath.c)
  endif(HAVE_FLEXPATH)

  if(HAVE_ICEE)
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} read/read_icee.c)
  endif()

  if(HAVE_NSSI)
    set(dist_libadiosreadf_nompi_a_SOURCES nssi/adios_nssi_args.x)
    set(nodist_libadiosreadf_nompi_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
    set(libadiosreadf_nompi_a_SOURCES ${libadiosreadf_nompi_a_SOURCES} ${dist_libadiosreadf_nompi_a_SOURCES} ${nodist_libadiosreadf_nompi_a_SOURCES})
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} nssi/adios_nssi_config.c read/read_nssi.c)
    set(libadiosreadf_nompi_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NAL")
    set(libadiosreadf_nompi_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")
    set(libadiosreadf_nompi_a_CFLAGS "${ADIOSLIB_CFLAGS}")
  endif(HAVE_NSSI)

  if(HAVE_NSSI)
    set(dist_libadiosreadf_nompi_a_SOURCES nssi/adios_nssi_args.x)
    set(nodist_libadiosreadf_nompi_a_SOURCES adios_nssi_args.c adios_nssi_args.h)
    set(libadiosreadf_nompi_a_SOURCES ${libadiosreadf_nompi_a_SOURCES} ${dist_libadiosreadf_nompi_a_SOURCES} ${nodist_libadiosreadf_nompi_a_SOURCES})
    set(FortranReadSeqLibSource ${FortranReadSeqLibSource} nssi/adios_nssi_config.c read/read_nssi.c)
    set(libadiosreadf_nompi_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS} ${MACRODEFFLAG}PTL_IFACE_CLIENT=CRAY_USER_NAL ${MACRODEFFLAG}PTL_IFACE_SERVER=CRAY_USER_NAL")
    set(libadiosreadf_nompi_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS} ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")
    set(libadiosreadf_nompi_a_CFLAGS "${ADIOSLIB_CFLAGS}")
  endif(HAVE_NSSI)

  set(libadiosreadf_nompi_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS} ${MACRODEFFLAG}_NOMPI ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_CPPFLAGS}")
  set(libadiosreadf_nompi_a_SOURCES ${libadiosreadf_nompi_a_SOURCES} ${FortranReadSeqLibSource} core/adiosf_read.c)
  set(libadiosreadf_nompi_v1_a_SOURCES ${FortranReadSeqLibSource} core/adiosf_read_v1.c)
  set(libadiosreadf_nompi_v1_a_CPPFLAGS "${libadiosreadf_nompi_a_CPPFLAGS}")
  set(libadiosreadf_nompi_v1_a_CFLAGS "${libadiosreadf_nompi_a_CFLAGS}")

  #add library adiosreadf_nompi and adiosreadf_nompi_v1
  add_library(adiosreadf_nompi ${libadiosreadf_nompi_a_SOURCES})
  set_target_properties(adiosreadf_nompi  PROPERTIES COMPILE_FLAGS "${libadiosreadf_nompi_a_CPPFLAGS} ${libadiosreadf_nompi_a_CFLAGS}")

  add_library(adiosreadf_nompi_v1 ${libadiosreadf_nompi_v1_a_SOURCES})
  set_target_properties(adiosreadf_nompi_v1 PROPERTIES COMPILE_FLAGS "${libadiosreadf_nompi_v1_a_CPPFLAGS} ${libadiosreadf_nompi_v1_a_CFLAGS}")

endif(BUILD_FORTRAN)


#start libadios_internal_nompi.a
set(libadios_internal_nompi_a_SOURCES core/mpidummy.c 
  core/adios_bp_v1.c 
  core/adios_endianness.c 
  core/bp_utils.c 
  core/adios_internals.c 
  ${transforms_common_SOURCES} 
  ${transforms_write_SOURCES} 
  ${query_C_SOURCES}
  core/buffer.c 
  core/adios_error.c 
  core/adios_logger.c 
  core/adios_timing.c 
  core/util.c core/strutil.c 
  core/qhashtbl.c 
  core/futils.c 
  core/adios_transport_hooks.c)

if(BUILD_WRITE)
  set(libadios_internal_nompi_a_SOURCES ${libadios_internal_nompi_a_SOURCES}
    core/adios_internals_mxml.c)
endif(BUILD_WRITE)

set(libadios_internal_nompi_a_CPPFLAGS "${MACRODEFFLAG}_INTERNAL ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSLIB_INT_CPPFLAGS} ${ADIOSLIB_CPPFLAGS}")

add_library(adios_internal_nompi ${libadios_internal_nompi_a_SOURCES})
set_target_properties(adios_internal_nompi PROPERTIES COMPILE_FLAGS "${libadios_internal_nompi_a_CPPFLAGS}")

add_subdirectory(core)
add_subdirectory(write)

install(FILES ${PROJECT_BINARY_DIR}/src/libadiosread.a
  ${PROJECT_BINARY_DIR}/src/libadiosread_nompi.a
  DESTINATION ${libdir})
if(BUILD_FORTRAN)
  install(FILES ${PROJECT_BINARY_DIR}/src/libadiosreadf.a
    ${PROJECT_BINARY_DIR}/src/libadiosreadf_nompi.a
    ${PROJECT_BINARY_DIR}/src/libadiosreadf_nompi_v1.a
    ${PROJECT_BINARY_DIR}/src/libadiosreadf_v1.a
    DESTINATION ${libdir})
endif()

if(BUILD_WRITE)
  install(FILES ${PROJECT_BINARY_DIR}/src/libadios.a
    ${PROJECT_BINARY_DIR}/src/libadios_nompi.a
    DESTINATION ${libdir})
  if(BUILD_FORTRAN)
    install(FILES ${PROJECT_BINARY_DIR}/src/libadiosf.a
      ${PROJECT_BINARY_DIR}/src/libadiosf_v1.a
      ${PROJECT_BINARY_DIR}/src/libadiosf_nompi.a
      ${PROJECT_BINARY_DIR}/src/libadiosf_nompi_v1.a
      DESTINATION ${libdir})
  endif()
endif()


message("setting sirius link/include paths")
add_library(splitdoubles STATIC IMPORTED)
set_target_properties(splitdoubles PROPERTIES
  IMPORTED_LOCATION "${SIRIUS_LIB}"
  IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
  )
target_link_libraries(adios splitdoubles)
target_link_libraries(adiosf splitdoubles)
target_link_libraries(adios_nompi splitdoubles)
target_link_libraries(adiosf_nompi splitdoubles)
#link_directories(${SIRIUS_LIB})
include_directories(${SIRIUS_INC})

#Note: ${MACRODEFFLAG}NOMPI chokes IBM's bgxlf compiler but it can pass ${MACRODEFFLAG}_NOMPI.

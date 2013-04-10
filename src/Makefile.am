INCLUDES=-I$(top_srcdir) 

AUTOMAKE_OPTIONS = no-dependencies
METASOURCES=AUTO

bin_PROGRAMS=
if HAVE_NSSI

bin_PROGRAMS+=nssi-staging-server

nssi_staging_server_SOURCES = nssi/nssi_staging_server.cpp nssi/aggregation.cpp
nssi_staging_server_CFLAGS=$(ADIOSLIB_CFLAGS)
nssi_staging_server_CPPFLAGS=$(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_CPPFLAGS) $(MACRODEFFLAG)PTL_IFACE_CLIENT=CRAY_USER_NAL $(MACRODEFFLAG)PTL_IFACE_SERVER=CRAY_USER_NAL
nssi_staging_server_LDFLAGS = $(ADIOSLIB_LDFLAGS) 
nssi_staging_server_LDADD = $(top_builddir)/src/libadios.a
nssi_staging_server_LDADD += $(ADIOSLIB_LDADD) $(NSSI_SERVER_LIBS)

bin_PROGRAMS+=nssi-coupling-server

nssi_coupling_server_SOURCES = nssi/nssi_coupling_server.cpp
nssi_coupling_server_CFLAGS=$(ADIOSLIB_CFLAGS)
nssi_coupling_server_CPPFLAGS=$(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_CPPFLAGS) $(MACRODEFFLAG)PTL_IFACE_CLIENT=CRAY_USER_NAL $(MACRODEFFLAG)PTL_IFACE_SERVER=CRAY_USER_NAL # $(MACRODEFFLAG)USE_TIMERS
nssi_coupling_server_LDFLAGS = $(ADIOSLIB_LDFLAGS) 
nssi_coupling_server_LDADD = $(top_builddir)/src/libadios.a
nssi_coupling_server_LDADD += $(ADIOSLIB_LDADD) $(NSSI_SERVER_LIBS)

endif # HAVE_NSSI

CLEANFILES= 
include_HEADERS=
nodist_include_HEADERS=
lib_LIBRARIES=

libadiosread_a_CPPFLAGS=
libadiosreadf_a_CPPFLAGS=
libadiosread_nompi_a_CPPFLAGS=
libadiosreadf_nompi_a_CPPFLAGS=

libadiosread_a_CFLAGS= 
libadiosreadf_a_CFLAGS= 
libadiosread_nompi_a_CFLAGS= 
libadiosreadf_nompi_a_CFLAGS= 

if BUILD_WRITE
lib_LIBRARIES += libadios.a libadios_nompi.a
libadios_a_SOURCES = core/adios.c \
                     core/common_adios.c \
                     core/adios_internals.c \
                     core/adios_internals_mxml.c \
                     core/buffer.c \
                     core/adios_bp_v1.c  \
                     core/adios_endianness.c \
                     core/bp_utils.c \
                     core/futils.c \
                     core/adios_error.c \
                     core/adios_read.c \
                     core/adios_read_v1.c \
                     core/common_read.c \
                     core/globals.c \
                     core/adios_timing.c \
                     core/adios_read_hooks.c \
                     core/adios_transport_hooks.c \
                     core/adios_socket.c \
                     core/adios_logger.c \
                     core/util.c \
                     read/read_bp.c \
                     read/read_bp_staged.c \
                     read/read_bp_staged1.c \
                     write/adios_mpi.c \
                     write/adios_mpi_lustre.c \
                     write/adios_mpi_amr.c \
                     write/adios_posix.c \
                     write/adios_posix1.c 


libadios_nompi_a_SOURCES = core/adios.c \
                     core/common_adios.c \
                     core/adios_internals.c \
                     core/adios_internals_mxml.c \
                     core/buffer.c \
                     core/adios_bp_v1.c  \
                     core/adios_endianness.c \
                     core/bp_utils.c \
                     core/futils.c \
                     core/adios_error.c \
                     core/adios_read.c \
                     core/adios_read_v1.c \
                     core/common_read.c \
                     core/globals.c \
                     core/mpidummy.c \
                     core/adios_timing.c \
                     core/adios_read_hooks.c \
                     core/adios_transport_hooks.c \
                     core/adios_socket.c \
                     core/adios_logger.c \
                     core/util.c \
                     read/read_bp.c \
                     read/read_bp_staged.c \
                     read/read_bp_staged1.c \
                     write/adios_posix.c \
                     write/adios_posix1.c 



if BUILD_FORTRAN
lib_LIBRARIES += libadiosf.a libadiosf_v1.a
FortranLibSources = core/adiosf.c \
                     core/common_adios.c \
                     core/adios_internals.c \
                     core/adios_internals_mxml.c \
                     core/buffer.c \
                     core/adios_bp_v1.c  \
                     core/adios_endianness.c\
                     core/futils.c \
                     core/adios_error.c \
                     core/bp_utils.c \
                     core/common_read.c \
                     core/globals.c \
                     core/adios_timing.c \
                     core/adios_read_hooks.c \
                     core/adios_transport_hooks.c \
                     core/adios_socket.c \
                     core/adios_logger.c \
                     core/util.c \
                     read/read_bp.c \
                     read/read_bp_staged.c \
                     read/read_bp_staged1.c \
                     write/adios_mpi.c \
                     write/adios_mpi_lustre.c \
                     write/adios_mpi_amr.c \
                     write/adios_posix.c \
                     write/adios_posix1.c 


endif # BUILD_FORTRAN

libadios_a_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_CPPFLAGS) # $(MACRODEFFLAG)USE_TIMERS
libadios_a_CFLAGS = $(ADIOSLIB_CFLAGS) 
libadios_nompi_a_CPPFLAGS = $(MACRODEFFLAG)_NOMPI $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_SEQ_CPPFLAGS) # $(MACRODEFFLAG)USE_TIMERS
libadios_nompi_a_CFLAGS = $(MACRODEFFLAG)_NOMPI $(ADIOSLIB_SEQ_CFLAGS) 
libadiosf_a_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_CPPFLAGS) # $(MACRODEFFLAG)USE_TIMERS
libadiosf_a_CFLAGS = $(ADIOSLIB_CFLAGS) 

if HAVE_DATASPACES
libadios_a_SOURCES += write/adios_dataspaces.c read/read_dataspaces.c
libadios_nompi_a_SOURCES += write/adios_dataspaces.c read/read_dataspaces.c
if BUILD_FORTRAN
FortranLibSources += write/adios_dataspaces.c read/read_dataspaces.c
endif
endif

if HAVE_DIMES
libadios_a_SOURCES += write/adios_dimes.c read/read_dimes.c
libadios_nompi_a_SOURCES += write/adios_dimes.c read/read_dimes.c
if BUILD_FORTRAN
FortranLibSources += write/adios_dimes.c read/read_dimes.c
endif
endif

if HAVE_DATATAP
libadios_a_SOURCES += write/adios_datatap.c read/read_datatap.c
libadios_nompi_a_SOURCES += write/adios_datatap.c read/read_datatap.c
if BUILD_FORTRAN
FortranLibSources += write/adios_datatap.c read/read_datatap.c
endif
endif 


if HAVE_PHDF5
libadios_a_SOURCES += write/adios_phdf5.c 
PHDF5_CPPFLAGS += $(MACRODEFFLAG)PHDF5 $(MACRODEFFLAG)H5_USE_16_API 
libadios_a_CPPFLAGS += $(PHDF5_CPPFLAGS)
if BUILD_FORTRAN
FortranLibSources += write/adios_phdf5.c 
libadiosf_a_CPPFLAGS += $(PHDF5_CPPFLAGS)
endif
endif

if HAVE_NC4PAR
libadios_a_SOURCES += write/adios_nc4.c 
if BUILD_FORTRAN
FortranLibSources += write/adios_nc4.c 
endif
endif

if HAVE_NSSI
dist_libadios_a_SOURCES = nssi/adios_nssi_args.x
nodist_libadios_a_SOURCES = adios_nssi_args.c adios_nssi_args.h
libadios_a_SOURCES += nssi/adios_nssi_config.c write/adios_nssi.c nssi/adios_nssi_filter.c read/read_nssi.c
libadios_a_CPPFLAGS += $(MACRODEFFLAG)PTL_IFACE_CLIENT=CRAY_USER_NAL $(MACRODEFFLAG)PTL_IFACE_SERVER=CRAY_USER_NAL
if BUILD_FORTRAN
nodist_libadiosf_a_SOURCES = adios_nssi_args.c adios_nssi_args.h
dist_libadiosf_a_SOURCES = nssi/adios_nssi_args.x
FortranLibSources += nssi/adios_nssi_config.c write/adios_nssi.c nssi/adios_nssi_filter.c read/read_nssi.c
libadiosf_a_CPPFLAGS += $(MACRODEFFLAG)PTL_IFACE_CLIENT=CRAY_USER_NAL $(MACRODEFFLAG)PTL_IFACE_SERVER=CRAY_USER_NAL
endif

# The next line forces automake to build adios_nssi_args.[ch] before compiling adios_nssi_config.c
# Note, the order of dist_*SOURCES and nodist_*SOURCES above should but does not solve this alone.
adios_nssi_config.c: adios_nssi_args.c

adios_nssi_args.h:  $(srcdir)/nssi/adios_nssi_args.x Makefile.am
	rm -f $@
#	@echo "#include <lwfs/lwfs_config.h>" >> $@
#	@echo "#ifndef LWFS_HAVE_STRUCT_SOCKADDR_IN" >> $@
#	@echo "#define LWFS_HAVE_STRUCT_SOCKADDR_IN" >> $@
#	@echo "struct sockaddr_in { " >> $@
#	@echo '  int a;' >> $@
#	@echo '};' >> $@
#	@echo '#endif' >> $@;
	rpcgen -h $(srcdir)/nssi/adios_nssi_args.x >> $@

adios_nssi_args.c: adios_nssi_args.h 
	rm -f adios_nssi_args.c
	rpcgen -c $(srcdir)/nssi/adios_nssi_args.x \
		| grep -v "include \"\.\." > adios_nssi_args.c

CLEANFILES += adios_nssi_args.c adios_nssi_args.h

endif # HAVE_NSSI


if RESEARCH_TRANSPORTS
libadios_a_SOURCES += write/adios_mpi_stagger.c \
                      write/adios_mpi_aggregate.c \
                      write/adios_provenance.c \
                      write/adios_mpi_cio.c \
                      write/adios_mpi_stripe.c \
                      write/adios_mpi_amr1.c \
                      write/adios_adaptive.c 
if BUILD_FORTRAN
FortranLibSources += write/adios_mpi_stagger.c \
                       write/adios_mpi_aggregate.c \
                       write/adios_provenance.c \
                       write/adios_mpi_cio.c \
                       write/adios_mpi_stripe.c \
                       write/adios_mpi_amr1.c \
                       write/adios_adaptive.c 
endif
else
libadios_a_CPPFLAGS += $(MACRODEFFLAG)NO_RESEARCH_TRANSPORTS
libadios_nompi_a_CPPFLAGS += $(MACRODEFFLAG)NO_RESEARCH_TRANSPORTS
libadiosf_a_CPPFLAGS += $(MACRODEFFLAG)NO_RESEARCH_TRANSPORTS
endif  # RESEARCH_TRANSPORTS

# Build two different Fortran library, for the two read APIs
if BUILD_FORTRAN
libadiosf_a_SOURCES = $(FortranLibSources) core/adiosf_read.c
libadiosf_v1_a_SOURCES = $(FortranLibSources) core/adiosf_read_v1.c
libadiosf_v1_a_CPPFLAGS = $(libadiosf_a_CPPFLAGS) 
libadiosf_v1_a_CFLAGS = $(libadiosf_a_CFLAGS) 

libadiosf_a_SOURCES += core/adiosf_defs_mod.f90 \
                       core/adiosf_write_mod.f90 \
                       core/adiosf_read_mod.f90 

libadiosf_v1_a_SOURCES += core/adiosf_defs_mod.f90 \
                          core/adiosf_write_mod.f90 

#nodist_libadiosreadf_a_SOURCES += adios_write.mod
nodist_include_HEADERS += adios_write_mod.mod 
CLEANFILES += adios_write_mod.mod

adios_write_mod.mod: core/adiosf_write_mod.f90 adios_defs_mod.mod
	$(FC) -c core/adiosf_write_mod.f90

endif


endif   # BUILD_WRITE

include_HEADERS += public/adios.h \
                   public/adios_types.h \
                   public/adios_read.h \
                   public/adios_error.h \
                   public/adios_mpi.h \
                   public/mpidummy.h \
                   public/adios_read_v1.h \
                   public/adios_read_v1_defs.h \
                   public/adios_read_v2.h \
                   public/adios_selection.h


lib_LIBRARIES += libadiosread.a
libadiosread_a_SOURCES = core/adios_bp_v1.c \
                      core/adios_endianness.c \
                      core/bp_utils.c \
                      core/futils.c \
                      core/adios_error.c \
                      core/adios_read.c \
                      core/adios_read_v1.c \
                      core/common_read.c \
                      core/globals.c \
                      core/adios_read_hooks.c \
                      core/adios_logger.c \
                      core/util.c \
                      read/read_bp.c \
                      read/read_bp_staged.c \
                      read/read_bp_staged1.c 
if HAVE_DATASPACES
libadiosread_a_SOURCES += read/read_dataspaces.c
endif
if HAVE_DIMES
libadiosread_a_SOURCES += read/read_dimes.c
endif
if HAVE_DATATAP
libadiosread_a_SOURCES += read/read_datatap.c
endif
if HAVE_NSSI
dist_libadiosread_a_SOURCES = nssi/adios_nssi_args.x
nodist_libadiosread_a_SOURCES = adios_nssi_args.c adios_nssi_args.h
libadiosread_a_SOURCES += nssi/adios_nssi_config.c read/read_nssi.c
libadiosread_a_CPPFLAGS += $(MACRODEFFLAG)PTL_IFACE_CLIENT=CRAY_USER_NAL $(MACRODEFFLAG)PTL_IFACE_SERVER=CRAY_USER_NAL
#libadiosread_a_CPPFLAGS += $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_CPPFLAGS) # $(MACRODEFFLAG)USE_TIMERS
libadiosread_a_CFLAGS   += $(ADIOSLIB_CFLAGS) 
endif

libadiosread_a_CPPFLAGS += $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_CPPFLAGS)

if BUILD_FORTRAN
lib_LIBRARIES += libadiosreadf.a libadiosreadf_v1.a
FortranReadLibSource = core/adios_bp_v1.c \
                      core/adios_endianness.c \
                      core/bp_utils.c \
                      core/futils.c \
                      core/adios_error.c \
                      core/common_read.c \
                      core/globals.c \
                      core/adios_read_hooks.c \
                      core/adios_logger.c \
                      core/util.c \
                      read/read_bp.c \
                      read/read_bp_staged.c \
                      read/read_bp_staged1.c 
                      
if HAVE_DATASPACES
FortranReadLibSource += read/read_dataspaces.c
endif
if HAVE_DIMES
FortranReadLibSource += read/read_dimes.c
endif
if HAVE_DATATAP
FortranReadLibSource += read/read_datatap.c
endif
if HAVE_NSSI
dist_libadiosreadf_a_SOURCES = nssi/adios_nssi_args.x
nodist_libadiosreadf_a_SOURCES = adios_nssi_args.c adios_nssi_args.h
FortranReadLibSource += nssi/adios_nssi_config.c read/read_nssi.c
libadiosreadf_a_CPPFLAGS += $(MACRODEFFLAG)PTL_IFACE_CLIENT=CRAY_USER_NAL $(MACRODEFFLAG)PTL_IFACE_SERVER=CRAY_USER_NAL
libadiosreadf_a_CFLAGS   += $(ADIOSLIB_CFLAGS) 
endif

libadiosreadf_a_CPPFLAGS += $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_CPPFLAGS)

libadiosreadf_a_SOURCES = $(FortranReadLibSource) core/adiosf_read.c
libadiosreadf_v1_a_SOURCES = $(FortranReadLibSource) core/adiosf_read_v1.c
libadiosreadf_v1_a_CPPFLAGS = $(libadiosreadf_a_CPPFLAGS)
libadiosreadf_v1_a_CFLAGS = $(libadiosreadf_a_CFLAGS)

libadiosreadf_a_SOURCES += core/adiosf_defs_mod.f90 \
                           core/adiosf_read_mod.f90

#nodist_libadiosreadf_a_SOURCES += adios_read_mod.mod
nodist_include_HEADERS += adios_read_mod.mod adios_defs_mod.mod
CLEANFILES += adios_read_mod.mod adios_defs_mod.mod

adios_read_mod.mod: core/adiosf_read_mod.f90 adios_defs_mod.mod 
	$(FC) -c core/adiosf_read_mod.f90

adios_defs_mod.mod: core/adiosf_defs_mod.f90
	$(FC) -c core/adiosf_defs_mod.f90

endif

lib_LIBRARIES += libadiosread_nompi.a
libadiosread_nompi_a_SOURCES = core/mpidummy.c\
                      core/adios_bp_v1.c \
                      core/adios_endianness.c \
                      core/bp_utils.c \
                      core/futils.c \
                      core/adios_error.c \
                      core/adios_read.c \
                      core/adios_read_v1.c \
                      core/common_read.c \
                      core/adios_logger.c \
                      core/buffer.c \
                      core/globals.c \
                      core/adios_read_hooks.c \
                      core/util.c \
                      read/read_bp.c 

if HAVE_DATASPACES
libadiosread_nompi_a_SOURCES += read/read_dataspaces.c
endif
if HAVE_DIMES
libadiosread_nompi_a_SOURCES += read/read_dimes.c
endif
if HAVE_DATATAP
libadiosread_nompi_a_SOURCES += read/read_datatap.c
endif
if HAVE_NSSI
dist_libadiosread_nompi_a_SOURCES = nssi/adios_nssi_args.x
nodist_libadiosread_nompi_a_SOURCES = adios_nssi_args.c adios_nssi_args.h
libadiosread_nompi_a_SOURCES += nssi/adios_nssi_config.c read/read_nssi.c
libadiosread_nompi_a_CPPFLAGS += $(MACRODEFFLAG)PTL_IFACE_CLIENT=CRAY_USER_NAL $(MACRODEFFLAG)PTL_IFACE_SERVER=CRAY_USER_NAL
libadiosread_nompi_a_CFLAGS   += $(ADIOSLIB_CFLAGS) 
endif
libadiosread_nompi_a_CPPFLAGS  += $(MACRODEFFLAG)_NOMPI $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_CPPFLAGS)

# Note: utils/adios_lint needs to link with adios_internals_mxml.o which depends
#       on adios_internals.c, buffer.c and adios_transport_hooks.c. 
#       We put these objects in libadiosread_nompi
#       so we do not need to make another _nompi lib just for adios_lint
if BUILD_FORTRAN
lib_LIBRARIES += libadiosreadf_nompi.a libadiosreadf_nompi_v1.a
FortranReadSeqLibSource = core/mpidummy.c\
                          core/adios_bp_v1.c \
                          core/adios_endianness.c \
                          core/bp_utils.c \
                          core/futils.c \
                          core/adios_error.c \
                          core/adios_logger.c \
                          core/common_read.c \
                          core/globals.c \
                          core/adios_read_hooks.c \
                          core/util.c \
                          read/read_bp.c 
if HAVE_DATASPACES
FortranReadSeqLibSource += read/read_dataspaces.c
endif
if HAVE_DIMES
FortranReadSeqLibSource += read/read_dimes.c
endif
if HAVE_DATATAP
FortranReadSeqLibSource += read/read_datatap.c
endif
if HAVE_NSSI
dist_libadiosreadf_nompi_a_SOURCES = nssi/adios_nssi_args.x
nodist_libadiosreadf_nompi_a_SOURCES = adios_nssi_args.c adios_nssi_args.h
FortranReadSeqLibSource += nssi/adios_nssi_config.c read/read_nssi.c
libadiosreadf_nompi_a_CPPFLAGS += $(MACRODEFFLAG)PTL_IFACE_CLIENT=CRAY_USER_NAL $(MACRODEFFLAG)PTL_IFACE_SERVER=CRAY_USER_NAL
libadiosreadf_nompi_a_CPPFLAGS += $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_CPPFLAGS) # $(MACRODEFFLAG)USE_TIMERS
libadiosreadf_nompi_a_CFLAGS   += $(ADIOSLIB_CFLAGS) 
endif

libadiosreadf_nompi_a_CPPFLAGS += $(MACRODEFFLAG)_NOMPI $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_CPPFLAGS)
libadiosreadf_nompi_a_SOURCES = $(FortranReadSeqLibSource) core/adiosf_read.c
libadiosreadf_nompi_v1_a_SOURCES = $(FortranReadSeqLibSource) core/adiosf_read_v1.c
libadiosreadf_nompi_v1_a_CPPFLAGS = $(libadiosreadf_nompi_a_CPPFLAGS)
libadiosreadf_nompi_v1_a_CFLAGS = $(libadiosreadf_nompi_a_CFLAGS)

endif  #BUILD_FORTRAN

noinst_LIBRARIES = libadios_internal_nompi.a
libadios_internal_nompi_a_SOURCES = core/mpidummy.c \
                                    core/adios_bp_v1.c \
                                    core/adios_endianness.c \
                                    core/bp_utils.c \
                                    core/adios_internals.c \
                                    core/buffer.c \
                                    core/adios_error.c \
                                    core/adios_logger.c \
                                    core/adios_internals_mxml.c \
                                    core/adios_timing.c \
                                    core/util.c \
                                    core/futils.c \
                                    core/adios_transport_hooks.c 


libadios_internal_nompi_a_CPPFLAGS = $(MACRODEFFLAG)_INTERNAL $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_CPPFLAGS)
#Note: $(MACRODEFFLAG)NOMPI chokes IBM's bgxlf compiler but it can pass $(MACRODEFFLAG)_NOMPI. 

if HAVE_DMALLOC
libadiosread_a_CPPFLAGS += $(MACRODEFFLAG)DMALLOC 
libadiosread_nompi_a_CPPFLAGS += $(MACRODEFFLAG)DMALLOC 
#if BUILD_FORTRAN
#libadiosreadf_a_CPPFLAGS += $(MACRODEFFLAG)DMALLOC $(DMALLOC_CPPFLAGS)
#libadiosreadf_nompi_a_CPPFLAGS += $(MACRODEFFLAG)DMALLOC $(DMALLOC_CPPFLAGS)
#endif
endif

CC=$(MPICC)
CXX=$(MPICXX)


#bin_PROGRAMS+=adios_flags
#nodist_adios_flags_SOURCES = adios_flags.h 
#adios_flags_SOURCES = adios_flags.c 
#adios_flags.h:  
#	rm -f $@
#	@echo "#define ADIOS_CFLAGS \"$(libadios_a_CPPFLAGS)\"" >> $@
#	@echo "#define ADIOS_LDFLAGS \"$(libadios_a_LDFLAGS)\"" >> $@
#	@echo "#define ADIOSF_CFLAGS \"$(libadiosf_a_CPPFLAGS)\"" >> $@
#	@echo "#define ADIOSF_LDFLAGS \"$(libadiosf_a_LDFLAGS)\"" >> $@
#	@echo "#define ADIOSREAD_CFLAGS \"$(libadiosread_a_CPPFLAGS)\"" >> $@
#	@echo "#define ADIOSREAD_LDFLAGS \"$(libadiosread_a_LDFLAGS)\"" >> $@
#	@echo "#define ADIOSREADF_CFLAGS \"$(libadiosreadf_a_CPPFLAGS)\"" >> $@
#	@echo "#define ADIOSREADF_LDFLAGS \"$(libadiosreadf_a_LDFLAGS)\"" >> $@
#	@echo "#define ADIOSREAD_NOMPI_CFLAGS \"$(libadiosread_nompi_a_CPPFLAGS)\"" >> $@
#	@echo "#define ADIOSREAD_NOMPI_LDFLAGS \"$(libadiosread_nompi_a_LDFLAGS)\"" >> $@
#	@echo "#define ADIOSREADF_NOMPI_CFLAGS \"$(libadiosreadf_nompi_a_CPPFLAGS)\"" >> $@
#	@echo "#define ADIOSREADF_NOMPI_LDFLAGS \"$(libadiosreadf_nompi_a_LDFLAGS)\"" >> $@
#adios_flags.c: adios_flags.h


EXTRA_DIST = core/adios_bp_v1.h core/adios_endianness.h \
             core/adios_internals.h core/adios_internals_mxml.h core/adios_logger.h \
             core/adios_read_hooks.h core/adios_socket.h core/adios_timing.h \
             core/adios_transport_hooks.h \
             core/bp_types.h core/bp_utils.h core/buffer.h core/common_adios.h \
             core/common_read.h core/futils.h core/globals.h core/ds_metadata.h \
             core/util.h \
             nssi/adios_nssi_config.h nssi/aggregation.h nssi/io_timer.h


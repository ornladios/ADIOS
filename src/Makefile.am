INCLUDES=-I$(top_srcdir) 

AUTOMAKE_OPTIONS = no-dependencies
METASOURCES=AUTO

bin_PROGRAMS=
if HAVE_NSSI
bin_PROGRAMS+=nssi-staging-server

nssi_staging_server_SOURCES = nssi_staging_server.cpp aggregation.cpp
nssi_staging_server_LDADD = $(top_builddir)/src/libadios.a
nssi_staging_server_LDFLAGS = $(MXML_LDFLAGS) 
nssi_staging_server_LDADD += $(MXML_LIBS) 
nssi_staging_server_CFLAGS=
nssi_staging_server_CPPFLAGS=

nssi_staging_server_CFLAGS += $(NSSI_CFLAGS)
nssi_staging_server_CPPFLAGS += $(NSSI_CPPFLAGS) -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL
nssi_staging_server_LDFLAGS += $(NSSI_LDFLAGS)
nssi_staging_server_LDADD += $(NSSI_LIBS) $(NSSI_SERVER_LIBS)

if HAVE_NC4PAR
nssi_staging_server_LDFLAGS += $(NC4PAR_LDFLAGS)
nssi_staging_server_LDADD +=  $(NC4PAR_LIBS)
endif
if HAVE_PHDF5
nssi_staging_server_LDFLAGS += $(PHDF5_LDFLAGS)
nssi_staging_server_LDADD +=  $(PHDF5_LIBS)
endif
if HAVE_DART
nssi_staging_server_CPPFLAGS += $(DART_CPPFLAGS)
nssi_staging_server_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
nssi_staging_server_LDADD += $(DART_LIBS)
endif

endif # HAVE_NSSI

include_HEADERS=
lib_LIBRARIES=

if BUILD_WRITE
lib_LIBRARIES += libadios.a
libadios_a_SOURCES = adios.c \
                     common_adios.c \
                     adios_transport_hooks.c \
                     adios_mpi.c \
                     adios_mpi_stripe.c \
                     adios_mpi_stripe2.c \
                     adios_mpi_amr.c \
                     adios_mpi_amr1.c \
                     adios_posix.c \
                     adios_posix1.c \
                     adios_internals.c \
                     adios_internals_mxml.c \
                     buffer.c \
                     adios_socket.c \
                     adios_bp_v1.c  \
                     adios_endianness.c \
                     bp_utils.c \
                     futils.c \
                     adios_error.c \
                     adios_read.c \
                     adios_read_hooks.c \
                     common_read.c \
                     read_bp.c \
                     globals.c
if BUILD_FORTRAN
lib_LIBRARIES += libadiosf.a
libadiosf_a_SOURCES = adiosf.c \
                     common_adios.c \
                     adios_transport_hooks.c \
                     adios_mpi.c \
                     adios_mpi_stripe.c \
                     adios_mpi_stripe2.c \
                     adios_mpi_amr.c \
                     adios_mpi_amr1.c \
                     adios_posix.c \
                     adios_posix1.c \
                     adios_internals.c \
                     adios_internals_mxml.c \
                     buffer.c \
                     adios_socket.c \
                     adios_bp_v1.c  \
                     adios_endianness.c\
                     futils.c \
                     adios_error.c \
                     bp_utils.c \
                     adiosf_read.c \
                     adios_read_hooks.c \
                     common_read.c \
                     read_bp.c \
                     globals.c
endif # BUILD_FORTRAN

libadios_a_CPPFLAGS = $(MXML_CPPFLAGS) # -DUSE_TIMERS
libadiosf_a_CPPFLAGS = $(MXML_CPPFLAGS)

if HAVE_DART
libadios_a_SOURCES += adios_dart.c read_dart.c
libadios_a_CPPFLAGS += $(DART_CPPFLAGS)
if BUILD_FORTRAN
libadiosf_a_SOURCES += adios_dart.c read_dart.c
libadiosf_a_CPPFLAGS += $(DART_CPPFLAGS)
endif
endif


if HAVE_DATATAP
libadios_a_SOURCES += adios_datatap.c
if BUILD_FORTRAN
libadiosf_a_SOURCES += adios_datatap.c
endif
endif 


if HAVE_PHDF5
libadios_a_SOURCES += adios_phdf5.c 
PHDF5_CPPFLAGS += -DPHDF5 -DH5_USE_16_API 
libadios_a_CPPFLAGS += $(PHDF5_CPPFLAGS)
if BUILD_FORTRAN
libadiosf_a_SOURCES += adios_phdf5.c 
libadiosf_a_CPPFLAGS += $(PHDF5_CPPFLAGS)
endif
endif

if HAVE_NC4PAR
libadios_a_SOURCES += adios_nc4.c 
libadios_a_CPPFLAGS += $(NC4PAR_CPPFLAGS)
if BUILD_FORTRAN
libadiosf_a_SOURCES += adios_nc4.c 
libadiosf_a_CPPFLAGS += $(NC4PAR_CPPFLAGS)
endif
endif

if HAVE_NSSI
libadios_a_SOURCES += adios_nssi_config.c adios_nssi.c adios_nssi_filter.c
dist_libadios_a_SOURCES = adios_nssi_args.x
nodist_libadios_a_SOURCES = adios_nssi_args.c adios_nssi_args.h
libadios_a_CPPFLAGS += $(NSSI_CPPFLAGS) -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL
if BUILD_FORTRAN
libadiosf_a_SOURCES += adios_nssi_config.c adios_nssi.c adios_nssi_filter.c
nodist_libadiosf_a_SOURCES = adios_nssi_args.c adios_nssi_args.h
dist_libadiosf_a_SOURCES = adios_nssi_args.x
libadiosf_a_CPPFLAGS += $(NSSI_CPPFLAGS) -DPTL_IFACE_CLIENT=CRAY_USER_NAL -DPTL_IFACE_SERVER=CRAY_USER_NAL
endif

adios_nssi_args.h:  $(srcdir)/adios_nssi_args.x Makefile.am
	rm -f $@
#	@echo "#include <lwfs/lwfs_config.h>" >> $@
#	@echo "#ifndef LWFS_HAVE_STRUCT_SOCKADDR_IN" >> $@
#	@echo "#define LWFS_HAVE_STRUCT_SOCKADDR_IN" >> $@
#	@echo "struct sockaddr_in { " >> $@
#	@echo '  int a;' >> $@
#	@echo '};' >> $@
#	@echo '#endif' >> $@;
	rpcgen -h $(srcdir)/adios_nssi_args.x >> $@

adios_nssi_args.c: adios_nssi_args.h 
	rm -f adios_nssi_args.c
	rpcgen -c $(srcdir)/adios_nssi_args.x \
		| grep -v "include \"\.\." > adios_nssi_args.c
endif # HAVE_NSSI


if RESEARCH_TRANSPORTS
libadios_a_SOURCES += adios_mpi_stagger.c \
                      adios_mpi_aggregate.c \
                      adios_provenance.c \
                      adios_mpi_cio.c \
                      adios_adaptive.c 
if BUILD_FORTRAN
libadiosf_a_SOURCES += adios_mpi_stagger.c \
                       adios_mpi_aggregate.c \
                       adios_provenance.c \
                       adios_mpi_cio.c \
                       adios_adaptive.c 
endif
else
libadios_a_CPPFLAGS += -DNO_RESEARCH_TRANSPORTS
libadiosf_a_CPPFLAGS += -DNO_RESEARCH_TRANSPORTS
endif  # RESEARCH_TRANSPORTS

endif   # BUILD_WRITE

include_HEADERS += adios.h \
                   adios_types.h \
                   adios_read.h \
                   adios_read_hooks.h \
                   mpidummy.h \
                   io_timer.h

lib_LIBRARIES += libadiosread.a
libadiosread_a_SOURCES = adios_bp_v1.c \
                      adios_endianness.c \
                      bp_utils.c \
                      futils.c \
                      adios_error.c \
                      adios_read.c \
                      adios_read_hooks.c \
                      common_read.c \
                      read_bp.c \
                      globals.c
if HAVE_DART
libadiosread_a_SOURCES += read_dart.c
endif
libadiosread_a_CPPFLAGS = $(DART_CPPFLAGS)

if BUILD_FORTRAN
lib_LIBRARIES += libadiosreadf.a
libadiosreadf_a_SOURCES = adios_bp_v1.c \
                      adios_endianness.c \
                      bp_utils.c \
                      futils.c \
                      adios_error.c \
                      adiosf_read.c \
                      adios_read_hooks.c \
                      common_read.c \
                      read_bp.c \
                      globals.c
                      
if HAVE_DART
libadiosreadf_a_SOURCES += read_dart.c
endif
libadiosreadf_a_CPPFLAGS = $(DART_CPPFLAGS)

# adiosf_read.mod is not really usable, 
# so it is not used and installed right now
#libadiosread_a_SOURCES += adiosf_read_mod.F90
#include_HEADERS += adiosf_read.mod 
endif

lib_LIBRARIES += libadiosread_nompi.a
libadiosread_nompi_a_SOURCES = mpidummy.c\
                      adios_bp_v1.c \
                      adios_endianness.c \
                      bp_utils.c \
                      futils.c \
                      adios_error.c \
                      adios_read.c \
                      adios_read_hooks.c \
                      common_read.c \
                      read_bp.c \
                      adios_internals.c \
                      buffer.c \
                      adios_transport_hooks.c \
                      globals.c

if HAVE_DART
libadiosread_nompi_a_SOURCES += read_dart.c
endif
# Note: utils/adios_lint needs to link with adios_internals_mxml.o which depends
#       on adios_internals.c, buffer.c and adios_transport_hooks.c. 
#       We put these objects in libadiosread_nompi
#       so we do not need to make another _nompi lib just for adios_lint
if BUILD_FORTRAN
lib_LIBRARIES += libadiosreadf_nompi.a
libadiosreadf_nompi_a_SOURCES = mpidummy.c\
                                adios_bp_v1.c \
                                adios_endianness.c \
                                bp_utils.c \
                                futils.c \
                                adios_error.c \
                                adiosf_read.c \
                                adios_read_hooks.c \
                                common_read.c \
                                read_bp.c 
if HAVE_DART
libadiosreadf_nompi_a_SOURCES += read_dart.c
endif
endif

noinst_LIBRARIES = libadios_internal_nompi.a
libadios_internal_nompi_a_SOURCES = mpidummy.c \
                                    adios_bp_v1.c \
                                    adios_endianness.c \
                                    bp_utils.c \
                                    adios_internals.c \
                                    buffer.c \
                                    adios_transport_hooks.c \
                                    adios_error.c \
                                    adios_internals_mxml.c


libadiosread_nompi_a_CPPFLAGS = -D_NOMPI $(DART_CPPFLAGS)
libadiosreadf_nompi_a_CPPFLAGS = -D_NOMPI $(DART_CPPFLAGS)
libadios_internal_nompi_a_CPPFLAGS = -D_NOMPI $(MXML_CPPFLAGS)
#Note: -DNOMPI chokes IBM's bgxlf compiler but it can pass -D_NOMPI. 

CC=$(MPICC)
CXX=$(MPICXX)

if HAVE_DMALLOC
libadiosread_a_CPPFLAGS += -DDMALLOC $(DMALLOC_CPPFLAGS)
libadiosread_nompi_a_CPPFLAGS += -DDMALLOC $(DMALLOC_CPPFLAGS)
#if BUILD_FORTRAN
#libadiosreadf_a_CPPFLAGS += -DDMALLOC $(DMALLOC_CPPFLAGS)
#libadiosreadf_nompi_a_CPPFLAGS += -DDMALLOC $(DMALLOC_CPPFLAGS)
#endif
endif

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


EXTRA_DIST = adios_bp_v1.h adios_bp_v2.h adios_endianness.h adios_error.h \
             adios_internals.h adios_internals_v2.h adios_internals_mxml.h \
             adios_nssi_config.h adios_socket.h adios_transport_hooks.h aggregation.h \
             bp_types.h bp_utils.h buffer.h common_adios.h common_read.h futils.h globals.h 


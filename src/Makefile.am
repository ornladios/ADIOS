INCLUDES=-I$(top_srcdir) 

AUTOMAKE_OPTIONS = no-dependencies
METASOURCES=AUTO

include_HEADERS=
lib_LIBRARIES=

if BUILD_WRITE
lib_LIBRARIES += libadios.a
libadios_a_SOURCES = adios.c \
                     common_adios.c \
                     adios_transport_hooks.c \
                     adios_mpi.c \
                     adios_mpi_stripe2.c \
                     adios_mpi_cio.c \
                     adios_mpi_stripe.c \
                     adios_mpi_stagger.c \
                     adios_mpi_aggregate.c \
                     adios_posix.c \
                     adios_posix_ascii.c \
                     adios_vtk.c \
                     adios_provenance.c \
                     adios_internals.c \
                     adios_socket.c \
                     adios_bp_v1.c  \
                     adios_bp_v2.c  \
                     adios_endianness.c\
                     adios_adaptive.c\
                     adios_datatap.c
#  		      adios_dart.c
if BUILD_FORTRAN
libadios_a_SOURCES += adiosf.c \
                      futils.c \
                      adios_error.c
endif

if HAVE_PHDF5
lib_LIBRARIES += libadios_phdf5.a
libadios_phdf5_a_SOURCES = adios_phdf5.c
PHDF5_CPPFLAGS += -DPHDF5 -DH5_USE_16_API 
libadios_phdf5_a_CPPFLAGS = $(MXML_CPPFLAGS)  $(PHDF5_CPPFLAGS)
else
libadios_a_SOURCES += adios_phdf5.c 
endif

#if DATATAP
#libadios_a_SOURCES += adios_datatap.c
#endif

libadios_a_CPPFLAGS = $(MXML_CPPFLAGS)  $(PHDF5_CPPFLAGS)

endif   # BUILD_WRITE

include_HEADERS += adios.h \
                   adios_types.h

lib_LIBRARIES += libadios_int.a
libadios_int_a_SOURCES= adios_internals.c \
                        adios_transport_hooks_dummy.c \
                        adios_bp_v1.c \
                        adios_endianness.c 
libadios_int_a_CPPFLAGS = $(MXML_CPPFLAGS)
#libadios_int_a_LIBADD = $(MXML_LIBS)

lib_LIBRARIES += libadios_int_nompi.a
libadios_int_nompi_a_SOURCES=adios_internals_nompi.c \
                        adios_transport_hooks_dummy.c \
                        adios_bp_v1.c \
                        adios_endianness.c 
libadios_int_nompi_a_CPPFLAGS = $(MXML_CPPFLAGS)
#libadios_int_nompi_a_LIBADD = $(MXML_LIBS)

lib_LIBRARIES += libadiosread.a
libadiosread_a_SOURCES = adios_bp_v1.c \
		      adios_endianness.c \
		      bp_utils.c \
                      futils.c \
                      adios_error.c \
		      adios_read.c
#		      adios_internals.c 
if BUILD_FORTRAN
libadiosread_a_SOURCES += adiosf_read.c adiosf_read_mod.F90
include_HEADERS += adiosf_read.mod 
#noinst_LIBRARIES = libadiosf_read_donotuse.a
#libadiosf_read_donotuse_a_SOURCES = adiosf_read.F90
endif
libadiosread_a_CPPFLAGS = $(MXML_CPPFLAGS)
#libadiosread_a_LIBADD = 

lib_LIBRARIES += libadiosread_nompi.a
libadiosread_nompi_a_SOURCES = mpidummy.c\
                      adios_bp_v1.c \
		      adios_endianness.c \
		      bp_utils.c \
                      futils.c \
                      adios_error.c \
		      adios_read.c 
if BUILD_FORTRAN
libadiosread_nompi_a_SOURCES += adiosf_read.c adiosf_read_mod.F90
endif

libadiosread_nompi_a_CPPFLAGS = $(MXML_CPPFLAGS) -DNOMPI 
include_HEADERS += adios_internals.h \
                   adios_transport_hooks.h 
include_HEADERS += adios_read.h \
                   mpidummy.h
CC=$(MPICC)

if HAVE_DMALLOC
libadiosread_nompi_a_CPPFLAGS += -DDMALLOC $(DMALLOC_CPPFLAGS)
libadiosread_a_CPPFLAGS += -DDMALLOC $(DMALLOC_CPPFLAGS)
endif

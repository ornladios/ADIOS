INCLUDES=-I$(top_srcdir) 

AUTOMAKE_OPTIONS = no-dependencies
METASOURCES=AUTO

include_HEADERS=
lib_LIBRARIES=

lib_LIBRARIES += libadios.a
libadios_a_SOURCES = adios.c \
                     adios_transport_hooks.c \
                     adios_mpi.c \
                     adios_mpi_stripe.c \
                     adios_mpi_cio.c \
                     adios_mpi_nwu.c \
                     adios_mpi_stagger.c \
                     adios_mpi_aggregate.c \
		     adios_posix.c \
                     adios_posix_ascii.c \
                     adios_vtk.c \
		     adios_provenance.c \
		     adios_internals.c \
		     adios_socket.c \
                     adios_bp_v1.c  \
                     adios_endianness.c\
		     adios_datatap.c
#		     adios_dart.c
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

lib_LIBRARIES += libbpread.a
libbpread_a_SOURCES = adios_bp_v1.c \
		      adios_endianness.c \
		      bp_utils.c \
		      bp_read.c
#		      adios_internals.c \
libbpread_a_CPPFLAGS = $(MXML_CPPFLAGS)
#libbpread_a_LIBADD = 

include_HEADERS += adios_internals.h \
                   adios_transport_hooks.h
include_HEADERS += bp_read.h
CC=$(MPICC)

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
                     adios_internals_mxml.c \
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
                   adios_types.h \
                   adios_read.h \
                   mpidummy.h

lib_LIBRARIES += libadiosread.a
libadiosread_a_SOURCES = adios_bp_v1.c \
		      adios_endianness.c \
		      bp_utils.c \
                      futils.c \
                      adios_error.c \
		      adios_read.c \
		      adios_internals.c 
if BUILD_FORTRAN
libadiosread_a_SOURCES += adiosf_read.c adiosf_read_mod.F90
include_HEADERS += adiosf_read.mod 
endif

lib_LIBRARIES += libadiosread_nompi.a
libadiosread_nompi_a_SOURCES = mpidummy.c\
                      adios_bp_v1.c \
		      adios_endianness.c \
		      bp_utils.c \
                      futils.c \
                      adios_error.c \
		      adios_read.c \
		      adios_internals.c \
                      adios_transport_hooks.c   
# Note: utils/adios_lint needs to link with adios_internals_mxml.o which depends
#       on adios_transport_hooks.c. We put this latter one in libadiosread_nompi
#       so we do not need to make another _nompi lib just for adios_lint
if BUILD_FORTRAN
libadiosread_nompi_a_SOURCES += adiosf_read.c adiosf_read_mod.F90
endif

libadiosread_nompi_a_CPPFLAGS = -DNOMPI
CC=$(MPICC)

if HAVE_DMALLOC
libadiosread_a_CPPFLAGS = -DDMALLOC $(DMALLOC_CPPFLAGS)
libadiosread_nompi_a_CPPFLAGS += -DDMALLOC $(DMALLOC_CPPFLAGS)
endif

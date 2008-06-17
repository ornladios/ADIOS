INCLUDES=-I$(top_srcdir)

AUTOMAKE_OPTIONS = no-dependencies
METASOURCES=AUTO

include_HEADERS=
lib_LIBRARIES=

lib_LIBRARIES += libadios.a
libadios_a_SOURCES = adios.c \
                     adios_transport_hooks.c \
                     adios_mpi.c \
                     adios_posix.c \
                     adios_posix_ascii.c \
                     adios_vtk.c
libadios_a_CPPFLAGS = $(MXML_CPPFLAGS)
include_HEADERS += adios.h \
                   adios_types.h

lib_LIBRARIES += libadios_int.a
libadios_int_a_SOURCES = adios_internals.c \
                         bw-utils.c \
                         br-utils.c \
                         binpack-utils.c \
                         adios_transport_hooks_dummy.c
libadios_int_a_CPPFLAGS = $(MXML_CPPFLAGS)
#libadios_int_a_LIBADD = $(MXML_LIBS)
include_HEADERS += adios_internals.h \
                   adios_transport_hooks.h \
                   br-utils.h \
                   bw-utils.h \
                   binpack-utils.h \
                   binpack-general.h



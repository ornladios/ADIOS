INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ascii

bp2ascii_SOURCES = bp2ascii.c \
                 hw-utils.c
bp2ascii_LDFLAGS = $(MXML_LDFLAGS)
bp2ascii_CPPFLAGS = $(MXML_CPPFLAGS)
bp2ascii_LDFLAGS += -L/apps/HDF5/1.6.5/ewok_par/pgi625/lib
bp2ascii_CPPFLAGS += -I/apps/HDF5/1.6.5/ewok_par/pgi625/include

bp2ascii_LDADD = $(top_builddir)/src/libadios.a
bp2ascii_LDADD += $(top_builddir)/src/libadios_int.a
bp2ascii_LDADD += $(MXML_LIBS)
bp2ascii_LDADD += -lhdf5 -lz

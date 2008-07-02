INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ncd

bp2ncd_SOURCES = bp2ncd.c

bp2ncd_LDFLAGS = $(MXML_LDFLAGS)
bp2ncd_CPPFLAGS = $(MXML_CPPFLAGS)

# change these to use proper autoconf settings
bp2ncd_LDFLAGS += -L/apps/HDF5/1.6.5/ewok_par/pgi625/lib
bp2ncd_CPPFLAGS += -I/apps/HDF5/1.6.5/ewok_par/pgi625/include

bp2ncd_LDADD = $(top_builddir)/src/libadios.a
bp2ncd_LDADD += $(top_builddir)/src/libadios_int.a
bp2ncd_LDADD += $(MXML_LIBS)
bp2ncd_LDADD += -lhdf5 -lz

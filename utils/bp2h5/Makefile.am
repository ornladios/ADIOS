INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2h5

bp2h5_SOURCES = bp2h5.c \
                 hw-utils.c
#                 parse_config.c
bp2h5_LDFLAGS = $(MXML_LDFLAGS)
#bp2h5_CPPFLAGS = $(MXML_CPPFLAGS)
bp2h5_CPPFLAGS = -DH5_USE_16_API
bp2h5_CPPFLAGS += $(MXML_CPPFLAGS)
# change these to use proper autoconf settings
bp2h5_LDFLAGS += $(HDF5_LDFLAGS)
bp2h5_CPPFLAGS += $(HDF5_CPPFLAGS)

bp2h5_LDADD = $(top_builddir)/src/libadios_int.a
bp2h5_LDADD += $(MXML_LIBS)
bp2h5_LDADD += -lhdf5 -lz

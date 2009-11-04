INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2h5

bp2h5_SOURCES = bp2h5.c hw-utils.c

bp2h5_CPPFLAGS = -DH5_USE_16_API
# change these to use proper autoconf settings
bp2h5_LDFLAGS = $(HDF5_LDFLAGS)
bp2h5_CPPFLAGS += $(HDF5_CPPFLAGS)

bp2h5_LDADD = $(top_builddir)/src/libadiosread_nompi.a
bp2h5_LDADD += -lhdf5 -lz

if HAVE_DMALLOC
bp2h5_CPPFLAGS += -DDMALLOC  $(DMALLOC_CPPFLAGS)
bp2h5_LDFLAGS +=  $(DMALLOC_LDFLAGS)
bp2h5_LDADD += -ldmalloc
endif

EXTRA_DIST = hw-utils.h


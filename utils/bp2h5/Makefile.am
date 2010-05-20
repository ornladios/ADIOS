INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2h5

bp2h5_SOURCES = bp2h5.c hw-utils.c

bp2h5_CPPFLAGS = -DH5_USE_16_API
# change these to use proper autoconf settings
bp2h5_LDFLAGS = -D_NOMPI $(HDF5_LDFLAGS) -static-libtool-libs
bp2h5_CPPFLAGS += $(HDF5_CPPFLAGS)

bp2h5_LDADD = $(top_builddir)/src/libadios_internal_nompi.a

if HAVE_DMALLOC
bp2h5_CPPFLAGS += -DDMALLOC  $(DMALLOC_CPPFLAGS)
bp2h5_LDFLAGS +=  $(DMALLOC_LDFLAGS)
bp2h5_LDADD += -ldmalloc
endif

CC=$(MPICC)

EXTRA_DIST = hw-utils.h


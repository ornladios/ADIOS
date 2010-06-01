INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ncd

bp2ncd_SOURCES = bp2ncd.c

# change these to use proper autoconf settings
bp2ncd_LDFLAGS = $(NETCDF_LDFLAGS)
bp2ncd_CPPFLAGS = -D_NOMPI $(NETCDF_CPPFLAGS) $(PHDF5_CPPFLAGS)

bp2ncd_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bp2ncd_LDADD += $(NETCDF_LIBS) $(PHDF5_LIBS) $(HDF5_LDFLAGS)

if HAVE_DMALLOC
bp2ncd_CPPFLAGS += -DDMALLOC  $(DMALLOC_CPPFLAGS)
bp2ncd_LDFLAGS +=  $(DMALLOC_LDFLAGS)
bp2ncd_LDADD += -ldmalloc
endif

EXTRA_DIST = bp2ncd.h
CC=$(MPICC)
CXX=$(MPICXX)

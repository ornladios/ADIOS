INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = adios_timeindex

adios_timeindex_SOURCES = adios_timeindex.c
adios_timeindex_LDADD = $(top_builddir)/src/libadios.a
adios_timeindex_LDFLAGS = $(MXML_LDFLAGS)
adios_timeindex_LDADD += $(MXML_LIBS)
if HAVE_PHDF5
adios_timeindex_LDFLAGS += $(PHDF5_LDFLAGS)
adios_timeindex_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)

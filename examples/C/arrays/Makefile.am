INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = arrays_write arrays_read

arrays_write_SOURCES = arrays_write.c
arrays_write_LDADD = $(top_builddir)/src/libadios.a 
arrays_write_LDFLAGS = $(MXML_LDFLAGS)
arrays_write_LDADD += $(MXML_LIBS)
arrays_read_SOURCES = arrays_read.c
arrays_read_LDADD = $(top_builddir)/src/libadios.a 
arrays_read_LDFLAGS = $(MXML_LDFLAGS)
arrays_read_LDADD += $(MXML_LIBS)

if HAVE_PHDF5
arrays_write_LDFLAGS += $(PHDF5_LDFLAGS)
arrays_write_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
arrays_read_LDFLAGS += $(PHDF5_LDFLAGS)
arrays_read_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = config.xml config_array.xml job.pbs

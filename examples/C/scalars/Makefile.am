INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = scalars_write scalars_read

scalars_write_SOURCES = scalars_write.c
scalars_write_LDADD = $(top_builddir)/src/libadios.a 
scalars_write_LDFLAGS = $(MXML_LDFLAGS)
scalars_write_LDADD += $(MXML_LIBS)
scalars_read_SOURCES = scalars_read.c
scalars_read_LDADD = $(top_builddir)/src/libadios.a 
scalars_read_LDFLAGS = $(MXML_LDFLAGS)
scalars_read_LDADD += $(MXML_LIBS)

if HAVE_PHDF5
scalars_write_LDFLAGS += $(PHDF5_LDFLAGS)
scalars_write_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
scalars_read_LDFLAGS += $(PHDF5_LDFLAGS)
scalars_read_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = config.xml config_array.xml job.pbs

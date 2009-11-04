INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_types adios_types_array adios_read_vars

adios_types_SOURCES = adios_types.c
adios_types_LDADD = $(top_builddir)/src/libadios.a 
adios_types_LDFLAGS = $(MXML_LDFLAGS)
adios_types_LDADD += $(MXML_LIBS)
adios_types_array_SOURCES = adios_types_array.c
adios_types_array_LDADD = $(top_builddir)/src/libadios.a 
adios_types_array_LDFLAGS = $(MXML_LDFLAGS)
adios_types_array_LDADD += $(MXML_LIBS)
adios_read_vars_SOURCES = adios_read_vars.c
adios_read_vars_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_vars_LDFLAGS = $(MXML_LDFLAGS)
adios_read_vars_LDADD += $(MXML_LIBS)

if HAVE_PHDF5
adios_types_LDFLAGS += $(PHDF5_LDFLAGS)
adios_types_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
adios_types_array_LDFLAGS += $(PHDF5_LDFLAGS)
adios_types_array_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = config.xml config_array.xml job.pbs

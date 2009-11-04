INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_global adios_global_2files

adios_global_SOURCES = adios_global.c
adios_global_LDADD = $(top_builddir)/src/libadios.a 
adios_global_LDFLAGS = $(MXML_LDFLAGS)
adios_global_LDADD += $(MXML_LIBS)
adios_global_2files_SOURCES = adios_global_2files.c
adios_global_2files_LDADD = $(top_builddir)/src/libadios.a 
adios_global_2files_LDFLAGS = $(MXML_LDFLAGS)
adios_global_2files_LDADD += $(MXML_LIBS)
if HAVE_PHDF5
adios_global_LDFLAGS += $(PHDF5_LDFLAGS)
adios_global_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
adios_global_2files_LDFLAGS += $(PHDF5_LDFLAGS)
adios_global_2files_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_global.xml gread_temperature.ch gwrite_temperature.ch job.pbs


INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_globaltime adios_read_globaltime

adios_globaltime_SOURCES = adios_globaltime.c
adios_globaltime_LDADD = $(top_builddir)/src/libadios.a
adios_globaltime_LDFLAGS = $(MXML_LDFLAGS)
adios_globaltime_LDADD += $(MXML_LIBS)
adios_read_globaltime_SOURCES = adios_read_globaltime.c
adios_read_globaltime_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_globaltime_LDFLAGS = $(MXML_LDFLAGS)
adios_read_globaltime_LDADD += $(MXML_LIBS)

if HAVE_PHDF5
adios_globaltime_LDFLAGS += $(PHDF5_LDFLAGS)
adios_globaltime_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
if HAVE_DART
adios_globaltime_CPPFLAGS = $(DART_CPPFLAGS)
adios_globaltime_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_globaltime_LDADD += $(DART_LIBS)
adios_read_globaltime_CPPFLAGS = $(DART_CPPFLAGS)
adios_read_globaltime_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_read_globaltime_LDADD += $(DART_LIBS)
endif
                
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_globaltime.xml gread_temperature.ch gwrite_temperature.ch job.pbs

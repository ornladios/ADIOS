INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_global adios_global_2files adios_read_global

adios_global_SOURCES = adios_global.c
adios_global_LDADD = $(top_builddir)/src/libadios.a 
adios_global_LDFLAGS = $(MXML_LDFLAGS)
adios_global_LDADD += $(MXML_LIBS)
adios_global_2files_SOURCES = adios_global_2files.c
adios_global_2files_LDADD = $(top_builddir)/src/libadios.a 
adios_global_2files_LDFLAGS = $(MXML_LDFLAGS)
adios_global_2files_LDADD += $(MXML_LIBS)
adios_read_global_SOURCES = adios_read_global.c
adios_read_global_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_global_LDFLAGS = $(MXML_LDFLAGS)
adios_read_global_LDADD += $(MXML_LIBS)

if HAVE_NETCDF
adios_global_LDFLAGS += $(NETCDF_LDFLAGS)
adios_global_LDADD +=  $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
adios_global_2files_LDFLAGS += $(NETCDF_LDFLAGS)
adios_global_2files_LDADD +=  $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
endif
if HAVE_NSSI
adios_global_LDFLAGS += $(NSSI_LDFLAGS)
adios_global_LDADD +=  $(top_builddir)/src/libadios_nssi.a $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
adios_global_2files_LDFLAGS += $(NSSI_LDFLAGS)
adios_global_2files_LDADD +=  $(top_builddir)/src/libadios_nssi.a $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
endif
if HAVE_PHDF5
adios_global_LDFLAGS += $(PHDF5_LDFLAGS)
adios_global_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
adios_global_2files_LDFLAGS += $(PHDF5_LDFLAGS)
adios_global_2files_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
if HAVE_DART
adios_global_CPPFLAGS = $(DART_CPPFLAGS)
adios_global_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_global_LDADD += $(DART_LIBS)
adios_global_2files_CPPFLAGS = $(DART_CPPFLAGS)
adios_global_2files_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_global_2files_LDADD += $(DART_LIBS)
adios_read_global_CPPFLAGS = $(DART_CPPFLAGS)
adios_read_global_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_read_global_LDADD += $(DART_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_global.xml gread_temperature.ch gwrite_temperature.ch job.pbs


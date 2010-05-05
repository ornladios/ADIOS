INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_globaltime

adios_globaltime_SOURCES = adios_globaltime.F90
adios_globaltime_LDADD = $(top_builddir)/src/libadiosf.a 
adios_globaltime_LDFLAGS = $(MXML_LDFLAGS)
adios_globaltime_LDADD += $(MXML_LIBS)
if HAVE_NETCDF
adios_globaltime_LDFLAGS += $(NETCDF_LDFLAGS)
adios_globaltime_LDADD +=   $(NETCDF_LIBS)
endif
if HAVE_NSSI
adios_globaltime_LDFLAGS += $(NSSI_LDFLAGS)
adios_globaltime_LDADD +=   $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
endif
if HAVE_PHDF5
adios_globaltime_LDFLAGS += $(PHDF5_LDFLAGS)
adios_globaltime_LDADD +=   $(PHDF5_LIBS)
endif
if HAVE_DART
adios_globaltime_CPPFLAGS = $(DART_CPPFLAGS)
adios_globaltime_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_globaltime_LDADD += $(DART_LIBS)
endif

CLEANFILES = *.bp
CC=$(MPICC)
FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

EXTRA_DIST = adios_globaltime.xml gread_temperature.fh gwrite_temperature.fh job.pbs

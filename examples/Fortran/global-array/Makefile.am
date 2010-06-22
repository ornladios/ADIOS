INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_global adios_global_no_xml

adios_global_SOURCES = adios_global.F90
adios_global_LDADD = $(top_builddir)/src/libadiosf.a 
adios_global_LDFLAGS = $(MXML_LDFLAGS)
adios_global_LDADD += $(MXML_LIBS)

adios_global_no_xml_SOURCES = adios_global_no_xml.F90
adios_global_no_xml_LDADD = $(top_builddir)/src/libadiosf.a 
adios_global_no_xml_LDFLAGS = $(MXML_LDFLAGS)
adios_global_no_xml_LDADD += $(MXML_LIBS)

if HAVE_NC4PAR
adios_global_LDFLAGS += $(NC4PAR_LDFLAGS)
adios_global_LDADD +=   $(NC4PAR_LIBS)
adios_global_no_xml_LDFLAGS += $(NC4PAR_LDFLAGS)
adios_global_no_xml_LDADD +=   $(NC4PAR_LIBS)
endif
if HAVE_NSSI
adios_global_LDFLAGS += $(NSSI_LDFLAGS)
adios_global_LDADD +=   $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
adios_global_no_xml_LDFLAGS += $(NSSI_LDFLAGS)
adios_global_no_xml_LDADD +=   $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
endif

if HAVE_PHDF5
adios_global_LDFLAGS += $(PHDF5_LDFLAGS)
adios_global_LDADD +=   $(PHDF5_LIBS)
adios_global_no_xml_LDFLAGS += $(PHDF5_LDFLAGS)
adios_global_no_xml_LDADD +=   $(PHDF5_LIBS)
endif
if HAVE_DART
adios_global_CPPFLAGS = $(DART_CPPFLAGS)
adios_global_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_global_LDADD += $(DART_LIBS)
adios_global_no_xml_CPPFLAGS = $(DART_CPPFLAGS)
adios_global_no_xml_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_global_no_xml_LDADD += $(DART_LIBS)
endif

CLEANFILES = *.bp
CC=$(MPICC)
FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

EXTRA_DIST = adios_global.xml gread_temperature.fh gwrite_temperature.fh job.pbs

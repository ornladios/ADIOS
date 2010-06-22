INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = 1_nonadios_example 2_adios_write 3_adios_read 4_adios_nfiles

1_nonadios_example_SOURCES = 1_nonadios_example.c
1_nonadios_example_LDADD = $(top_builddir)/src/libadios.a 
1_nonadios_example_LDFLAGS = $(MXML_LDFLAGS)
1_nonadios_example_LDADD += $(MXML_LIBS)
2_adios_write_SOURCES = 2_adios_write.c
2_adios_write_LDADD = $(top_builddir)/src/libadios.a 
2_adios_write_LDFLAGS = $(MXML_LDFLAGS)
2_adios_write_LDADD += $(MXML_LIBS)
3_adios_read_SOURCES = 3_adios_read.c
3_adios_read_LDADD = $(top_builddir)/src/libadios.a 
3_adios_read_LDFLAGS = $(MXML_LDFLAGS)
3_adios_read_LDADD += $(MXML_LIBS)
4_adios_nfiles_SOURCES = 4_adios_nfiles.c
4_adios_nfiles_LDADD = $(top_builddir)/src/libadios.a 
4_adios_nfiles_LDFLAGS = $(MXML_LDFLAGS)
4_adios_nfiles_LDADD += $(MXML_LIBS)

if HAVE_NC4PAR
1_nonadios_example_LDFLAGS += $(NC4PAR_LDFLAGS)
1_nonadios_example_LDADD +=   $(NC4PAR_LIBS)
2_adios_write_LDFLAGS += $(NC4PAR_LDFLAGS)
2_adios_write_LDADD +=   $(NC4PAR_LIBS)
3_adios_read_LDFLAGS += $(NC4PAR_LDFLAGS)
3_adios_read_LDADD +=   $(NC4PAR_LIBS)
4_adios_nfiles_LDFLAGS += $(NC4PAR_LDFLAGS)
4_adios_nfiles_LDADD +=   $(NC4PAR_LIBS)
endif
if HAVE_NSSI
1_nonadios_example_LDFLAGS += $(NSSI_LDFLAGS)
1_nonadios_example_LDADD +=   $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
2_adios_write_LDFLAGS += $(NSSI_LDFLAGS)
2_adios_write_LDADD +=   $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
3_adios_read_LDFLAGS += $(NSSI_LDFLAGS)
3_adios_read_LDADD +=   $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
4_adios_nfiles_LDFLAGS += $(NSSI_LDFLAGS)
4_adios_nfiles_LDADD +=   $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
endif
if HAVE_PHDF5
1_nonadios_example_LDFLAGS += $(PHDF5_LDFLAGS)
1_nonadios_example_LDADD +=   $(PHDF5_LIBS)
2_adios_write_LDFLAGS += $(PHDF5_LDFLAGS)
2_adios_write_LDADD +=   $(PHDF5_LIBS)
3_adios_read_LDFLAGS += $(PHDF5_LDFLAGS)
3_adios_read_LDADD +=   $(PHDF5_LIBS)
4_adios_nfiles_LDFLAGS += $(PHDF5_LDFLAGS)
4_adios_nfiles_LDADD +=   $(PHDF5_LIBS)
endif
if HAVE_DART
1_nonadios_example_CPPFLAGS = $(DART_CPPFLAGS)
1_nonadios_example_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
1_nonadios_example_LDADD += $(DART_LIBS)
2_adios_write_CPPFLAGS = $(DART_CPPFLAGS)
2_adios_write_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
2_adios_write_LDADD += $(DART_LIBS)
3_adios_read_CPPFLAGS = $(DART_CPPFLAGS)
3_adios_read_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
3_adios_read_LDADD += $(DART_LIBS)
4_adios_nfiles_CPPFLAGS = $(DART_CPPFLAGS)
4_adios_nfiles_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
4_adios_nfiles_LDADD += $(DART_LIBS)
endif

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = config.xml gread_temperature.ch gwrite_temperature.ch

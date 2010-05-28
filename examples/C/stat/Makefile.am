INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = stat_write stat_read

stat_write_SOURCES = stat_write.c
stat_write_LDADD = $(top_builddir)/src/libadios.a 
stat_write_LDFLAGS = $(MXML_LDFLAGS)
stat_write_LDADD += $(MXML_LIBS)
stat_read_SOURCES = stat_read.c
stat_read_LDADD = $(top_builddir)/src/libadios.a 
stat_read_LDFLAGS = $(MXML_LDFLAGS)
stat_read_LDADD += $(MXML_LIBS)

if HAVE_NETCDF
stat_write_LDFLAGS += $(NETCDF_LDFLAGS)
stat_write_LDADD +=  $(NETCDF_LIBS)
stat_read_LDFLAGS += $(NETCDF_LDFLAGS)
stat_read_LDADD +=  $(NETCDF_LIBS)
endif
if HAVE_NSSI
stat_write_LDFLAGS += $(NSSI_LDFLAGS)
stat_write_LDADD +=  $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
stat_read_LDFLAGS += $(NSSI_LDFLAGS)
stat_read_LDADD +=  $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
endif
if HAVE_PHDF5
stat_write_LDFLAGS += $(PHDF5_LDFLAGS)
stat_write_LDADD +=  $(PHDF5_LIBS)
stat_read_LDFLAGS += $(PHDF5_LDFLAGS)
stat_read_LDADD +=  $(PHDF5_LIBS)
endif
if HAVE_DART
stat_write_CPPFLAGS = $(DART_CPPFLAGS)
stat_write_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
stat_write_LDADD += $(DART_LIBS)
stat_read_CPPFLAGS = $(DART_CPPFLAGS)
stat_read_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
stat_read_LDADD += $(DART_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = stat.xml gwrite_stat.ch gread_stat.ch job.pbs

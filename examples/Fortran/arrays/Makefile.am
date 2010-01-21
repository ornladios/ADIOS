INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = arrays_write arrays_read

arrays_write_SOURCES = arrays_write.F90
arrays_write_LDADD = $(top_builddir)/src/libadiosf.a 
arrays_write_LDFLAGS = $(MXML_LDFLAGS)
arrays_write_LDADD += $(MXML_LIBS)

arrays_read_SOURCES = arrays_read.F90
arrays_read_LDADD = $(top_builddir)/src/libadiosf.a 
arrays_read_LDFLAGS = $(MXML_LDFLAGS)
arrays_read_LDADD += $(MXML_LIBS)

if HAVE_NETCDF
arrays_write_LDFLAGS += $(NETCDF_LDFLAGS)
arrays_write_LDADD +=  $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
arrays_read_LDFLAGS += $(NETCDF_LDFLAGS)
arrays_read_LDADD +=  $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
endif
if HAVE_PHDF5
arrays_write_LDFLAGS += $(PHDF5_LDFLAGS)
arrays_write_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
arrays_read_LDFLAGS += $(PHDF5_LDFLAGS)
arrays_read_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
if HAVE_DART
arrays_write_CPPFLAGS = $(DART_CPPFLAGS)
arrays_write_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
arrays_write_LDADD += $(DART_LIBS)
arrays_read_CPPFLAGS = $(DART_CPPFLAGS)
arrays_read_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
arrays_read_LDADD += $(DART_LIBS)
endif

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = arrays.xml gwrite_arrays.fh gread_arrays.fh job.pbs

FC=$(MPIFC)

FCLINK = \
        $(LIBTOOL) --mode=link --tag F77 $(FC) \
        $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

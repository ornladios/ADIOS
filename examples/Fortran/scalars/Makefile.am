INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = scalars_write scalars_read

scalars_write_SOURCES = scalars_write.F90
scalars_write_LDADD = $(top_builddir)/src/libadiosf.a 
scalars_write_LDFLAGS = $(MXML_LDFLAGS)
scalars_write_LDADD += $(MXML_LIBS)

scalars_read_SOURCES = scalars_read.F90
scalars_read_LDADD = $(top_builddir)/src/libadiosf.a 
scalars_read_LDFLAGS = $(MXML_LDFLAGS)
scalars_read_LDADD += $(MXML_LIBS)

if HAVE_NETCDF
scalars_write_LDFLAGS += $(NETCDF_LDFLAGS)
scalars_write_LDADD +=  $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
scalars_read_LDFLAGS += $(NETCDF_LDFLAGS)
scalars_read_LDADD +=  $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
endif
if HAVE_PHDF5
scalars_write_LDFLAGS += $(PHDF5_LDFLAGS)
scalars_write_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
scalars_read_LDFLAGS += $(PHDF5_LDFLAGS)
scalars_read_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
if HAVE_DART
scalars_write_CPPFLAGS = $(DART_CPPFLAGS)
scalars_write_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
scalars_write_LDADD += $(DART_LIBS)
scalars_read_CPPFLAGS = $(DART_CPPFLAGS)
scalars_read_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
scalars_read_LDADD += $(DART_LIBS)
endif

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = scalars.xml gwrite_scalars.fh gread_scalars.fh job.pbs

FC=$(MPIFC)

FCLINK = \
        $(LIBTOOL) --mode=link --tag F77 $(FC) \
        $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

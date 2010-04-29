INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS=genarray genarray2D copyarray copyarray2D

genarray_SOURCES = genarray.F90
genarray_LDFLAGS = $(MXML_LDFLAGS)
genarray_LDADD = $(top_builddir)/src/libadiosf.a $(MXML_LIBS)
if HAVE_NETCDF
genarray_LDADD += $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
genarray_LDFLAGS += $(NETCDF_LDFLAGS)
endif
if HAVE_NSSI
genarray_LDADD += $(top_builddir)/src/libadios_nssi.a $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
genarray_LDFLAGS += $(NSSI_LDFLAGS)
endif
if HAVE_PHDF5
genarray_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
genarray_LDFLAGS += $(PHDF5_LDFLAGS)
endif
if HAVE_DART
genarray_CPPFLAGS = $(DART_CPPFLAGS)
genarray_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
genarray_LDADD += $(DART_LIBS)
endif

genarray2D_SOURCES = genarray2D.F90
genarray2D_LDFLAGS = $(MXML_LDFLAGS)
genarray2D_LDADD = $(top_builddir)/src/libadiosf.a $(MXML_LIBS)
if HAVE_NETCDF
genarray2D_LDADD += $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
genarray2D_LDFLAGS += $(NETCDF_LDFLAGS)
endif
if HAVE_NSSI
genarray2D_LDADD += $(top_builddir)/src/libadios_nssi.a $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
genarray2D_LDFLAGS += $(NSSI_LDFLAGS)
endif
if HAVE_PHDF5
genarray2D_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
genarray2D_LDFLAGS += $(PHDF5_LDFLAGS)
endif
if HAVE_DART
genarray2D_CPPFLAGS = $(DART_CPPFLAGS)
genarray2D_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
genarray2D_LDADD += $(DART_LIBS)
endif

copyarray_SOURCES = copyarray.F90
copyarray_LDFLAGS = $(MXML_LDFLAGS)
copyarray_LDADD = $(top_builddir)/src/libadiosf.a $(MXML_LIBS)
if HAVE_NETCDF
copyarray_LDADD += $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
copyarray_LDFLAGS += $(NETCDF_LDFLAGS)
endif
if HAVE_NSSI
copyarray_LDADD += $(top_builddir)/src/libadios_nssi.a $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
copyarray_LDFLAGS += $(NSSI_LDFLAGS)
endif
if HAVE_PHDF5
copyarray_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
copyarray_LDFLAGS += $(PHDF5_LDFLAGS)
endif
if HAVE_DART
copyarray_CPPFLAGS = $(DART_CPPFLAGS)
copyarray_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
copyarray_LDADD += $(DART_LIBS)
endif

copyarray2D_SOURCES = copyarray2D.F90
copyarray2D_LDFLAGS = $(MXML_LDFLAGS)
copyarray2D_LDADD = $(top_builddir)/src/libadiosf.a $(MXML_LIBS)
if HAVE_NETCDF
copyarray2D_LDADD += $(top_builddir)/src/libadios_nc4.a $(NETCDF_LIBS)
copyarray2D_LDFLAGS += $(NETCDF_LDFLAGS)
endif
if HAVE_NSSI
copyarray2D_LDADD += $(top_builddir)/src/libadios_nssi.a $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
copyarray2D_LDFLAGS += $(NSSI_LDFLAGS)
endif
if HAVE_PHDF5
copyarray2D_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
copyarray2D_LDFLAGS += $(PHDF5_LDFLAGS)
endif
if HAVE_DART
copyarray2D_CPPFLAGS = $(DART_CPPFLAGS)
copyarray2D_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
copyarray2D_LDADD += $(DART_LIBS)
endif

FC=$(MPIFC)

genarray.F90: gwrite_genarray.fh
copyarray.F90: gwrite_genarray.fh gread_genarray.fh

gwrite_genarray.fh: genarray3d.xml
	rm -f gwrite_genarray.fh gread_genarray.fh
	$(top_builddir)/utils/gpp/gpp.py genarray3d.xml
gread_genarray.fh: genarray3d.xml
	rm -f gwrite_genarray.fh gread_genarray.fh
	$(top_builddir)/utils/gpp/gpp.py genarray3d.xml

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

CLEANFILES = genarray_comm.mod genarray2d_comm.mod gwrite_genarray.fh gread_genarray.fh

EXTRA_DIST = genarray.xml genarray3d.xml job.ewok job2D.ewok

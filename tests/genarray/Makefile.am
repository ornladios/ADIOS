INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS=genarray genarray2D

genarray_SOURCES = genarray.F90
genarray_LDFLAGS = $(MXML_LDFLAGS)
genarray_LDADD = $(top_builddir)/src/libadiosf.a $(MXML_LIBS)
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
if HAVE_PHDF5
genarray2D_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
genarray2D_LDFLAGS += $(PHDF5_LDFLAGS)
endif
if HAVE_DART
genarray2D_CPPFLAGS = $(DART_CPPFLAGS)
genarray2D_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
genarray2D_LDADD += $(DART_LIBS)
endif

FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

CLEANFILES = genarray_comm.mod genarray2d_comm.mod

EXTRA_DIST = genarray.xml job.ewok job2D.ewok

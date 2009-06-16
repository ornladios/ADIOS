INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=genarray genarray2D

genarray_SOURCES = genarray.F90
genarray_LDFLAGS = $(MXML_LDFLAGS)
genarray_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS)
if HAVE_PHDF5
genarray_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
genarray_LDFLAGS += $(PHDF5_LDFLAGS)
endif

genarray2D_SOURCES = genarray2D.F90
genarray2D_LDFLAGS = $(MXML_LDFLAGS)
genarray2D_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS)
if HAVE_PHDF5
genarray2D_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
genarray2D_LDFLAGS += $(PHDF5_LDFLAGS)
endif

FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

CLEANFILES = genarray_comm.mod genarray2D_comm.mod


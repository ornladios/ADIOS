INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS=adios_test_f

#adios_test_f_SOURCES = adios_test_f.F90 adios_module.F90
adios_test_f_SOURCES = adios_test_f.F90
adios_test_f_LDFLAGS = $(MXML_LDFLAGS)
adios_test_f_LDADD = $(top_builddir)/src/libadiosf.a $(MXML_LIBS)
if HAVE_PHDF5
adios_test_f_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
adios_test_f_LDFLAGS += $(PHDF5_LDFLAGS)
endif

FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

EXTRA_DIST = config_fortran.xml

INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS=adios_test_f

#adios_test_f_SOURCES = adios_test_f.F90 adios_module.F90
adios_test_f_SOURCES = adios_test_f.F90
adios_test_f_LDFLAGS = $(MXML_LDFLAGS)
adios_test_f_LDADD = $(top_builddir)/src/libadiosf.a $(MXML_LIBS)
if HAVE_NETCDF
adios_test_f_LDADD += $(NETCDF_LIBS)
adios_test_f_LDFLAGS += $(NETCDF_LDFLAGS)
endif
if HAVE_NSSI
adios_test_f_LDADD += $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
adios_test_f_LDFLAGS += $(NSSI_LDFLAGS)
endif
if HAVE_PHDF5
adios_test_f_LDADD +=  $(PHDF5_LIBS)
adios_test_f_LDFLAGS += $(PHDF5_LDFLAGS)
endif
if HAVE_DART
adios_test_f_CPPFLAGS = $(DART_CPPFLAGS)
adios_test_f_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_test_f_LDADD += $(DART_LIBS)
endif

FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

EXTRA_DIST = config_fortran.xml

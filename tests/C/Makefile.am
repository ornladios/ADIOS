INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS=adios_test_c

adios_test_c_SOURCES = adios_test_c.c
adios_test_c_LDADD = $(top_builddir)/src/libadios.a
adios_test_c_LDFLAGS = $(MXML_LDFLAGS) 
adios_test_c_LDADD += $(MXML_LIBS) 
if HAVE_NETCDF
adios_test_c_LDFLAGS += $(NETCDF_LDFLAGS)
adios_test_c_LDADD +=  $(NETCDF_LIBS)
endif
if HAVE_NSSI
adios_test_c_LDADD += $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
adios_test_c_LDFLAGS += $(NSSI_LDFLAGS)
endif
if HAVE_PHDF5
adios_test_c_LDFLAGS += $(PHDF5_LDFLAGS)
adios_test_c_LDADD +=  $(PHDF5_LIBS)
endif
if HAVE_DART
adios_test_c_CPPFLAGS = $(DART_CPPFLAGS)
adios_test_c_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
adios_test_c_LDADD += $(DART_LIBS)
endif
CLEANFILES = restart.bp
CC=$(MPICC)

EXTRA_DIST = config_c.xml

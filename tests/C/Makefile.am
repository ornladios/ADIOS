INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS=adios_test_c

adios_test_c_SOURCES = adios_test_c.c
adios_test_c_CPPFLAGS = $(ADIOSLIB_CPPFLAGS)
adios_test_c_CFLAGS = $(ADIOSLIB_CFLAGS)
adios_test_c_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
adios_test_c_LDFLAGS = $(ADIOSLIB_LDFLAGS) 

CLEANFILES = restart.bp
CC=$(MPICC)

EXTRA_DIST = config_c.xml

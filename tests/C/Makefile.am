INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=adios_test_c

adios_test_c_SOURCES = adios_test_c.c
adios_test_c_LDFLAGS = $(MXML_LDFLAGS)
adios_test_c_LDADD = $(top_builddir)/src/libadios.a
adios_test_c_LDADD += $(MXML_LIBS)
CLEANFILES = restart.0

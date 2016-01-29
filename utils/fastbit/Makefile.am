AUTOMAKE_OPTIONS = no-dependencies

AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/core


bin_PROGRAMS = adios_index_fastbit

adios_index_fastbit_SOURCES = adios_index_fastbit.c 
adios_index_fastbit_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_CPPFLAGS) $(ADIOSLIB_CFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS)
adios_index_fastbit_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
adios_index_fastbit_LDFLAGS = $(ADIOSLIB_LDFLAGS)

CLEANFILES = *.idx


CC=$(MPICC)

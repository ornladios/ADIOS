INCLUDES = $(all_includes)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

#all-local:
#	test "$(srcdir)" = "$(builddir)" && $(MAKE) $(AM_MAKEFLAGS) all 
	 
noinst_PROGRAMS = adios_alac_query

adios_alac_query_SOURCES = adios_alac_query.c
adios_alac_query_LDADD = $(top_builddir)/src/libadios.a 
adios_alac_query_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_alac_query_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_alac_query.xml gread_temperature.ch gwrite_temperature.ch job.pbs


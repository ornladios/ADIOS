AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public -I$(top_srcdir)/src  -I$(top_srcdir)/src/public -I$(top_srcdir)/src/core

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=adios_lint

adios_lint_SOURCES = adios_lint.c
adios_lint_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
adios_lint_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS)
adios_lint_LDADD = $(top_builddir)/src/libadios_internal_nompi.a 
adios_lint_LDADD += $(ADIOSLIB_INT_LDADD)


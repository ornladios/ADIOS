AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src  -I$(top_builddir)/src/public  -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(top_srcdir)/src/core -I$(top_srcdir)/src/core/transforms

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpmeta

bpmeta_SOURCES = bpmeta.c
bpmeta_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
bpmeta_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS)
bpmeta_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bpmeta_LDADD += $(ADIOSLIB_INT_LDADD)


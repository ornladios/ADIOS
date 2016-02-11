AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public  -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(top_srcdir)/src/core

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = bp2ascii

bp2ascii_SOURCES = bp2ascii.c 
bp2ascii_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
bp2ascii_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS)
bp2ascii_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bp2ascii_LDADD += $(ADIOSLIB_INT_LDADD)


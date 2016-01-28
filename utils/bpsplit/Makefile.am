AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(top_srcdir)/src/core

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpsplit bpappend bpgettime

bpsplit_SOURCES = bpsplit.c
bpsplit_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
bpsplit_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS)
bpsplit_LDADD = $(top_builddir)/src/libadios_internal_nompi.a 
bpsplit_LDADD += $(ADIOSLIB_INT_LDADD)

bpappend_SOURCES = bpappend.c
bpappend_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
bpappend_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS)
bpappend_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bpappend_LDADD += $(ADIOSLIB_INT_LDADD)

bpgettime_SOURCES = bpgettime.c
bpgettime_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
bpgettime_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS)
bpgettime_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bpgettime_LDADD += $(ADIOSLIB_INT_LDADD)


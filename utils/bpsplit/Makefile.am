AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(top_srcdir)/src/core

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpsplit bpappend bpgettime

bpsplit_SOURCES = bpsplit.c
bpsplit_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bpsplit_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
bpsplit_LDADD = $(top_builddir)/src/libadios_nompi.a 
bpsplit_LDADD += $(ADIOSLIB_SEQ_LDADD)

bpappend_SOURCES = bpappend.c
bpappend_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bpappend_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
bpappend_LDADD = $(top_builddir)/src/libadios_nompi.a
bpappend_LDADD += $(ADIOSLIB_SEQ_LDADD)

bpgettime_SOURCES = bpgettime.c
bpgettime_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bpgettime_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
bpgettime_LDADD = $(top_builddir)/src/libadios_nompi.a
bpgettime_LDADD += $(ADIOSLIB_SEQ_LDADD)


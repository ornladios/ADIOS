AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(top_srcdir)/src/core

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ncd

bp2ncd_SOURCES = bp2ncd.c

if USE_PARALLEL_NETCDF_FOR_UTILS
    bp2ncd_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_CPPFLAGS)
    bp2ncd_CFLAGS = $(ADIOSLIB_CFLAGS) 
    bp2ncd_LDFLAGS = $(ADIOSLIB_LDFLAGS) -static-libtool-libs
    bp2ncd_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
else
    bp2ncd_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
    bp2ncd_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS) -static-libtool-libs
    bp2ncd_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
    bp2ncd_LDADD += $(ADIOSLIB_INT_LDADD)
endif

CC=$(BP2NCD_CC)

EXTRA_DIST = bp2ncd.h

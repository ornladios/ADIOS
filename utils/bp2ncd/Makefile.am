INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ncd

bp2ncd_SOURCES = bp2ncd.c

if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
    bp2ncd_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_CPPFLAGS)
    bp2ncd_CFLAGS = $(ADIOSLIB_CFLAGS) 
    bp2ncd_LDFLAGS = $(ADIOSLIB_LDFLAGS) -static-libtool-libs
    bp2ncd_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
else
    bp2ncd_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
    bp2ncd_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS) -static-libtool-libs
    bp2ncd_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
    bp2ncd_LDADD += $(ADIOSLIB_INT_LDADD)
endif

EXTRA_DIST = bp2ncd.h

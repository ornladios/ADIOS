INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpdump

bpdump_SOURCES = bpdump.c
bpdump_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
bpdump_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS)
bpdump_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bpdump_LDADD += $(ADIOSLIB_INT_LDADD)

if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
endif

INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpdump

bpdump_SOURCES = bpdump.c
bpdump_CPPFLAGS = $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bpdump_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS)
bpdump_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bpdump_LDADD += $(ADIOSLIB_SEQ_LDADD)

if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
endif

INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpsplit bpappend bpgettime

bpsplit_SOURCES = bpsplit.c
bpsplit_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bpsplit_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS)
bpsplit_LDADD = $(top_builddir)/src/libadios_internal_nompi.a 
bpsplit_LDADD += $(ADIOSLIB_SEQ_LDADD)

bpappend_SOURCES = bpappend.c
bpappend_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bpappend_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS)
bpappend_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bpappend_LDADD += $(ADIOSLIB_SEQ_LDADD)

bpgettime_SOURCES = bpgettime.c
bpgettime_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bpgettime_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS)
bpgettime_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bpgettime_LDADD += $(ADIOSLIB_SEQ_LDADD)

if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
endif


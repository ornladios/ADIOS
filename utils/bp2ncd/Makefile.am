INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ncd

bp2ncd_SOURCES = bp2ncd.c

bp2ncd_CPPFLAGS = $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bp2ncd_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS) -static-libtool-libs
bp2ncd_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bp2ncd_LDADD += $(ADIOSLIB_SEQ_LDADD)

EXTRA_DIST = bp2ncd.h
if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
endif

INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ascii

bp2ascii_SOURCES = bp2ascii.c 
bp2ascii_CPPFLAGS = $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bp2ascii_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS)
bp2ascii_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bp2ascii_LDADD += $(ADIOSLIB_SEQ_LDADD)

if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
endif

INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2bp

bp2bp_SOURCES = bp2bp.c
bp2bp_CPPFLAGS = $(ADIOSLIB_CPPFLAGS) $(ADIOSLIB_CFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) 
bp2bp_LDFLAGS = $(ADIOSLIB_LDFLAGS) 
bp2bp_LDADD =  $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)


#if USE_PARALLEL_COMPILER
CC=$(MPICC)
CXX=$(MPICXX)
#endif



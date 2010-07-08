INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=adios_lint

adios_lint_SOURCES = adios_lint.c
adios_lint_CPPFLAGS = $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
adios_lint_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS)
adios_lint_LDADD = $(top_builddir)/src/libadios_internal_nompi.a 
adios_lint_LDADD += $(ADIOSLIB_SEQ_LDADD)

if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
endif

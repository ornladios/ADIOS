INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=adios_lint

adios_lint_SOURCES = adios_lint.c
adios_lint_CPPFLAGS = $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_INT_CPPFLAGS) $(ADIOSLIB_INT_CFLAGS)
adios_lint_LDFLAGS = $(ADIOSLIB_INT_LDFLAGS)
adios_lint_LDADD = $(top_builddir)/src/libadios_internal_nompi.a 
adios_lint_LDADD += $(ADIOSLIB_INT_LDADD)

if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
endif

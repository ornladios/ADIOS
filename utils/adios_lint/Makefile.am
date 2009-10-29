INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=adios_lint

adios_lint_SOURCES = adios_lint.c
adios_lint_CPPFLAGS = $(MXML_CPPFLAGS)
adios_lint_LDFLAGS = $(MXML_LDFLAGS)
adios_lint_LDADD = $(top_builddir)/src/libadiosread_nompi.a $(top_builddir)/src/libadios_a-adios_internals_mxml.o
adios_lint_LDADD += $(MXML_LIBS)

if HAVE_DMALLOC
adios_lint_CPPFLAGS += -DDMALLOC  $(DMALLOC_CPPFLAGS)
adios_lint_LDFLAGS +=  $(DMALLOC_LDFLAGS)
adios_lint_LDADD += -ldmalloc
endif

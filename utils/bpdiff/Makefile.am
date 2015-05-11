AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpdiff

bpdiff_SOURCES = bpdiff.c decompose_block.c utils.c
bpdiff_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_CPPFLAGS) $(ADIOSREADLIB_CFLAGS)
bpdiff_LDFLAGS = $(ADIOSREADLIB_LDFLAGS) 
bpdiff_LDADD = $(top_builddir)/src/libadiosread.a 
bpdiff_LDADD += $(ADIOSREADLIB_LDADD)

CC=$(MPICC)
CXX=$(MPICXX)

EXTRA_DIST = decompose.h utils.h

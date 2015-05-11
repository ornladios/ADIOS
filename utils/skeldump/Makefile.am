AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(top_srcdir)/src/core

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = skeldump 

skeldump_SOURCES = skeldump.c 
skeldump_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_SEQ_CPPFLAGS) $(ADIOSREADLIB_SEQ_CFLAGS)
skeldump_LDFLAGS = $(ADIOSREADLIB_SEQ_LDFLAGS) 
skeldump_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
skeldump_LDADD += $(ADIOSREADLIB_SEQ_LDADD)

EXTRA_DIST = skeldump.h
#CC=$(MPICC)
#CXX=$(MPICXX)


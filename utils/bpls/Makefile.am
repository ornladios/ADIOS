AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpls

bpls_SOURCES = bpls.c 
bpls_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_SEQ_CPPFLAGS) $(ADIOSREADLIB_SEQ_CFLAGS)
bpls_LDFLAGS = $(ADIOSREADLIB_SEQ_LDFLAGS) 
bpls_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
bpls_LDADD += $(ADIOSREADLIB_SEQ_LDADD)

EXTRA_DIST = bpls.h
#CC=$(MPICC)
#CXX=$(MPICXX)


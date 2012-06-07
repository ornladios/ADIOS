INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I$(top_builddir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpls

bpls_SOURCES = bpls.c 
bpls_CPPFLAGS = -DADIOS_USE_READ_API_1 $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_SEQ_CPPFLAGS) $(ADIOSREADLIB_SEQ_CFLAGS)
bpls_LDFLAGS = $(ADIOSREADLIB_SEQ_LDFLAGS) 
bpls_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
bpls_LDADD += $(ADIOSREADLIB_SEQ_LDADD)

EXTRA_DIST = bpls.h
#CC=$(MPICC)
#CXX=$(MPICXX)


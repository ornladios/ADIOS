INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpls

bpls_SOURCES = bpls.c 
bpls_CFLAGS = -DNOMPI -g 
bpls_LDADD = $(top_builddir)/src/libadiosread_nompi.a 

if HAVE_DMALLOC
bpls_CPPFLAGS = -DDMALLOC  $(DMALLOC_CPPFLAGS)
bpls_LDFLAGS =  $(DMALLOC_LDFLAGS)
bpls_LDADD += -ldmalloc
endif

EXTRA_DIST = bpls.h


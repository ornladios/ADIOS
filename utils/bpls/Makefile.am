INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpls

bpls_SOURCES = bpls.c 
bpls_CFLAGS = -D_NOMPI -g 
bpls_LDFLAGS = 
bpls_LDADD = $(top_builddir)/src/libadiosread_nompi.a 

if HAVE_DMALLOC
bpls_CPPFLAGS = -DDMALLOC  $(DMALLOC_CPPFLAGS)
bpls_LDFLAGS +=  $(DMALLOC_LDFLAGS)
bpls_LDADD += -ldmalloc
endif

if HAVE_DART
bpls_CPPFLAGS = $(DART_CPPFLAGS)
bpls_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS) 
bpls_LDADD += $(DART_LIBS)
endif

EXTRA_DIST = bpls.h


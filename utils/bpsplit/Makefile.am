INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpsplit bpappend bpgettime

bpsplit_SOURCES = bpsplit.c
bpsplit_LDADD = $(top_builddir)/src/libadiosread_nompi.a

bpappend_SOURCES = bpappend.c
bpappend_LDADD = $(top_builddir)/src/libadiosread_nompi.a

bpgettime_SOURCES = bpgettime.c
bpgettime_LDADD = $(top_builddir)/src/libadiosread_nompi.a

if HAVE_DMALLOC
bpsplit_CPPFLAGS = -DDMALLOC  $(DMALLOC_CPPFLAGS)
bpappend_CPPFLAGS = -DDMALLOC  $(DMALLOC_CPPFLAGS)
bpgettime_CPPFLAGS = -DDMALLOC  $(DMALLOC_CPPFLAGS)
bpsplit_LDFLAGS =  $(DMALLOC_LDFLAGS)
bpappend_LDFLAGS =  $(DMALLOC_LDFLAGS)
bpgettime_LDFLAGS =  $(DMALLOC_LDFLAGS)
bpsplit_LDADD += -ldmalloc
bpappend_LDADD += -ldmalloc
bpgettime_LDADD += -ldmalloc
endif

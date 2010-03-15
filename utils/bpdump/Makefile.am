INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpdump bpdump_v2

bpdump_SOURCES = bpdump.c
bpdump_CPPFLAGS = -D_NOMPI
bpdump_LDADD = $(top_builddir)/src/libadios_internal_nompi.a

bpdump_v2_SOURCES = bpdump_v2.c
bpdump_v2_CPPFLAGS = -D_NOMPI
bpdump_v2_LDADD = $(top_builddir)/src/libadios_internal_nompi.a

if HAVE_DMALLOC
bpdump_CPPFLAGS += -DDMALLOC  $(DMALLOC_CPPFLAGS)
bpdump_LDFLAGS =  $(DMALLOC_LDFLAGS)
bpdump_LDADD += -ldmalloc
bpdump_v2_CPPFLAGS += -DDMALLOC  $(DMALLOC_CPPFLAGS)
bpdump_v2_LDFLAGS =  $(DMALLOC_LDFLAGS)
bpdump_v2_LDADD += -ldmalloc
endif

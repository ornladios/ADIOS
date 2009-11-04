INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpdump

bpdump_SOURCES = bpdump.c
bpdump_CPPFLAGS = -D_NOMPI
bpdump_LDADD = $(top_builddir)/src/libadiosread_nompi.a

if HAVE_DMALLOC
bpdump_CPPFLAGS += -DDMALLOC  $(DMALLOC_CPPFLAGS)
bpdump_LDFLAGS =  $(DMALLOC_LDFLAGS)
bpdump_LDADD += -ldmalloc
endif

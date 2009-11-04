INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ascii

bp2ascii_SOURCES = bp2ascii.c 

bp2ascii_LDADD = $(top_builddir)/src/libadiosread_nompi.a

if HAVE_DMALLOC
bp2ascii_CPPFLAGS = -DDMALLOC  $(DMALLOC_CPPFLAGS)
bp2ascii_LDFLAGS =  $(DMALLOC_LDFLAGS)
bp2ascii_LDADD += -ldmalloc
endif



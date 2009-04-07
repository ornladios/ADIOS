INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpdump

bpdump_SOURCES = bpdump.c
bpdump_LDFLAGS = $(MXML_LDFLAGS)
bpdump_LDADD = $(top_builddir)/src/libadios_int_nompi.a
bpdump_LDADD += $(MXML_LIBS)

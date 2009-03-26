INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpsplit bpappend

bpsplit_SOURCES = bpsplit.c
bpsplit_LDFLAGS = $(MXML_LDFLAGS)
bpsplit_LDADD = $(top_builddir)/src/libadios_int.a
bpsplit_LDADD += $(MXML_LIBS)

bpappend_SOURCES = bpappend.c
bpappend_LDFLAGS = $(MXML_LDFLAGS)
bpappend_LDADD = $(top_builddir)/src/libadios_int.a
bpappend_LDADD += $(MXML_LIBS)

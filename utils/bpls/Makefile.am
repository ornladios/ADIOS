INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bpls

bpls_SOURCES = bpls.c 
bpls_CFLAGS = -DNOMPI -g
bpls_LDFLAGS = $(MXML_LDFLAGS)
bpls_LDADD = $(top_builddir)/src/libadiosread_nompi.a
bpls_LDADD += $(MXML_LIBS)

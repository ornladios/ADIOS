INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = read_all

read_all_SOURCES = read_all.c
read_all_CPPFLAGS = -D_NOMPI
read_all_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
read_all_LDFLAGS = $(MXML_LDFLAGS)
read_all_LDADD += $(MXML_LIBS)
if HAVE_DART
read_all_CPPFLAGS += $(DART_CPPFLAGS)
read_all_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
read_all_LDADD += $(DART_LIBS)
endif

CC=$(MPICC)

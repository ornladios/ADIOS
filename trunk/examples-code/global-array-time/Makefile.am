INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = adios_globaltime

adios_globaltime_SOURCES = adios_globaltime.c
adios_globaltime_LDADD = $(top_builddir)/src/libadios.a
adios_globaltime_LDFLAGS = $(MXML_LDFLAGS)
adios_globaltime_LDADD += $(MXML_LIBS)
if HAVE_PHDF5
adios_globaltime_LDFLAGS += $(PHDF5_LDFLAGS)
adios_globaltime_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp

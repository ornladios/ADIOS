AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = adios_globaltime adios_globaltime_no_xml adios_read_globaltime adios_read_globaltime_streaming

adios_globaltime_SOURCES = adios_globaltime.c
adios_globaltime_LDADD = $(top_builddir)/src/libadios.a
adios_globaltime_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_globaltime_LDADD += $(ADIOSLIB_LDADD)
adios_globaltime_no_xml_SOURCES = adios_globaltime_no_xml.c
adios_globaltime_no_xml_LDADD = $(top_builddir)/src/libadios.a
adios_globaltime_no_xml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_globaltime_no_xml_LDADD += $(ADIOSLIB_LDADD)
adios_read_globaltime_SOURCES = adios_read_globaltime.c
adios_read_globaltime_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_globaltime_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
adios_read_globaltime_LDADD += $(ADIOSREADLIB_LDADD)
adios_read_globaltime_streaming_SOURCES = adios_read_globaltime_streaming.c
adios_read_globaltime_streaming_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_globaltime_streaming_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
adios_read_globaltime_streaming_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_globaltime.xml gread_restart.ch gwrite_restart.ch job.pbs

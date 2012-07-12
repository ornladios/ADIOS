INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I$(top_builddir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

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
adios_read_globaltime_CPPFLAGS = -DADIOS_USE_READ_API_1
adios_read_globaltime_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_globaltime_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_globaltime_LDADD += $(ADIOSLIB_LDADD)
adios_read_globaltime_streaming_SOURCES = adios_read_globaltime_streaming.c
adios_read_globaltime_streaming_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_globaltime_streaming_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_globaltime_streaming_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_globaltime.xml gread_restart.ch gwrite_restart.ch job.pbs

INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_globaltime adios_globaltime_no_xml adios_read_globaltime

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
adios_read_globaltime_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_globaltime_LDADD += $(ADIOSLIB_LDADD)
                
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_globaltime.xml gread_temperature.ch gwrite_temperature.ch job.pbs

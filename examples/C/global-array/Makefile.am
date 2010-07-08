INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_global adios_global_2files adios_read_global adios_global_no_xml adios_read_global_no_xml

adios_global_SOURCES = adios_global.c
adios_global_LDADD = $(top_builddir)/src/libadios.a 
adios_global_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_LDADD += $(ADIOSLIB_LDADD)
adios_global_2files_SOURCES = adios_global_2files.c
adios_global_2files_LDADD = $(top_builddir)/src/libadios.a 
adios_global_2files_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_2files_LDADD += $(ADIOSLIB_LDADD)
adios_read_global_SOURCES = adios_read_global.c
adios_read_global_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_global_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_global_LDADD += $(ADIOSLIB_LDADD)
adios_global_no_xml_SOURCES = adios_global_no_xml.c
adios_global_no_xml_LDADD = $(top_builddir)/src/libadios.a 
adios_global_no_xml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_no_xml_LDADD += $(ADIOSLIB_LDADD)
adios_read_global_no_xml_SOURCES = adios_read_global_no_xml.c
adios_read_global_no_xml_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_global_no_xml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_global_no_xml_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_global.xml gread_temperature.ch gwrite_temperature.ch job.pbs


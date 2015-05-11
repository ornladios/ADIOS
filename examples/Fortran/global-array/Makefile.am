AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = adios_global adios_global_no_xml no_xml_write_byid

adios_global_SOURCES = adios_global.F90
adios_global_LDADD = $(top_builddir)/src/libadiosf.a 
adios_global_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_LDADD += $(ADIOSLIB_LDADD)

adios_global_no_xml_SOURCES = adios_global_no_xml.F90
adios_global_no_xml_LDADD = $(top_builddir)/src/libadiosf.a 
adios_global_no_xml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_no_xml_LDADD += $(ADIOSLIB_LDADD)

no_xml_write_byid_SOURCES = no_xml_write_byid.F90
no_xml_write_byid_LDADD = $(top_builddir)/src/libadiosf.a
no_xml_write_byid_LDFLAGS = $(ADIOSLIB_LDFLAGS)
no_xml_write_byid_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

EXTRA_DIST = adios_global.xml gread_temperature.fh gwrite_temperature.fh job.pbs

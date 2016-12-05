AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = time_aggregation time_aggregation_no_xml

time_aggregation_SOURCES = time_aggregation.c
time_aggregation_LDADD = $(top_builddir)/src/libadios.a
time_aggregation_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
time_aggregation_LDADD += $(ADIOSLIB_LDADD)

time_aggregation_no_xml_SOURCES = time_aggregation_no_xml.c
time_aggregation_no_xml_LDADD = $(top_builddir)/src/libadios.a
time_aggregation_no_xml_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
time_aggregation_no_xml_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = time_aggregation.xml

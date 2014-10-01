INCLUDES = $(all_includes)
INCLUDES += $(ADIOSLIB_CPPFLAGS)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = build_indexed_dataset compute_expected_query_results
noinst_LIBRARIES = libparse_test_query_xml.a

build_indexed_dataset_SOURCES = build_indexed_dataset.c 
build_indexed_dataset_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD) 
build_indexed_dataset_LDFLAGS = $(ADIOSLIB_LDFLAGS)

libparse_test_query_xml_a_SOURCES = adios_query_xml_parse.c adios_query_xml_parse.h
libparse_test_query_xml_a_LIBADD = $(top_builddir)/src/libadios.a

compute_expected_query_results_SOURCES = compute_expected_query_results.c
compute_expected_query_results_LDADD = $(builddir)/libparse_test_query_xml.a $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD) 
compute_expected_query_results_LDFLAGS = $(ADIOSLIB_LDFLAGS)

CLEANFILES = *.bp
CC=$(MPICC)
#EXTRA_DIST = test_config_flex.xml 

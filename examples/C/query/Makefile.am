INCLUDES = $(all_includes)
INCLUDES += -I$(top_srcdir)/src/public -I$(top_srcdir)/src/query

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_read_ext_test

adios_read_ext_test_SOURCES = adios_read_ext_test.c
adios_read_ext_test_LDADD = $(top_builddir)/src/libadios.a 
adios_read_ext_test_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_ext_test_LDADD += $(ADIOSLIB_LDADD)

if HAVE_ALACRITY
noinst_PROGRAMS += adios_build_alac_index adios_read_alac_index adios_show_bp_file

adios_build_alac_index_SOURCES = adios_build_alac_index.c
adios_build_alac_index_LDADD = $(top_builddir)/src/libadios.a 
adios_build_alac_index_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_build_alac_index_LDADD += $(ADIOSLIB_LDADD)
adios_build_alac_index_CPPFLAGS = $(ALACRITY_CPPFLAGS)

adios_read_alac_index_SOURCES = adios_read_alac_index.c
adios_read_alac_index_LDADD = $(top_builddir)/src/libadios.a 
adios_read_alac_index_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_alac_index_LDADD += $(ADIOSLIB_LDADD)
adios_read_alac_index_CPPFLAGS = $(ALACRITY_CPPFLAGS)

adios_show_bp_file_SOURCES = adios_show_bp_file.c
adios_show_bp_file_LDADD = $(top_builddir)/src/libadios.a 
adios_show_bp_file_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_show_bp_file_LDADD += $(ADIOSLIB_LDADD)
endif # HAVE_ALACRITY

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_build_alac_index.xml gread_temperature.ch gwrite_temperature.ch job.pbs


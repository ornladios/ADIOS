AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public  -I$(top_srcdir)/src -I$(top_srcdir)/src/public $(ZFP_CPPFLAGS)

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = adios_global adios_global_2files adios_read_gpp adios_read_global test_xgc test_zfp test_f0_zfp adios_read_writeblock adios_read_chunk adios_global_no_xml adios_read_global_no_xml no_xml_write_byid read_no_xml_write_byid adios_global_aggregate_by_color

adios_global_SOURCES = adios_global.c
adios_global_LDADD = $(top_builddir)/src/libadios.a 
adios_global_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_LDADD += $(ADIOSLIB_LDADD)
adios_global_2files_SOURCES = adios_global_2files.c
adios_global_2files_LDADD = $(top_builddir)/src/libadios.a 
adios_global_2files_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_2files_LDADD += $(ADIOSLIB_LDADD)
adios_read_global_SOURCES = adios_read_global.c
adios_read_global_LDADD = $(top_builddir)/src/libadios.a 
adios_read_global_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_global_LDADD += $(ADIOSLIB_LDADD)
test_xgc_SOURCES = test_xgc.c
test_xgc_LDADD = $(top_builddir)/src/libadios.a
test_xgc_LDFLAGS = $(ADIOSLIB_LDFLAGS) 
test_xgc_LDADD += $(ADIOSLIB_LDADD)
test_zfp_SOURCES = test_zfp.c
test_zfp_LDADD = $(top_builddir)/src/libadios.a
test_zfp_LDFLAGS = $(ADIOSLIB_LDFLAGS)
test_zfp_LDADD += $(ADIOSLIB_LDADD)
test_f0_zfp_SOURCES = test_f0_zfp.c
test_f0_zfp_LDADD = $(top_builddir)/src/libadios.a
test_f0_zfp_LDFLAGS = $(ADIOSLIB_LDFLAGS)
test_f0_zfp_LDADD += $(ADIOSLIB_LDADD)

adios_read_writeblock_SOURCES = adios_read_writeblock.c
adios_read_writeblock_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_writeblock_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_writeblock_LDADD += $(ADIOSLIB_LDADD)
adios_read_chunk_SOURCES = adios_read_chunk.c
adios_read_chunk_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_chunk_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
adios_read_chunk_LDADD += $(ADIOSREADLIB_LDADD)
adios_read_gpp_SOURCES = adios_read_gpp.c
adios_read_gpp_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_gpp_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
adios_read_gpp_LDADD += $(ADIOSREADLIB_LDADD)
adios_global_no_xml_SOURCES = adios_global_no_xml.c
adios_global_no_xml_LDADD = $(top_builddir)/src/libadios.a 
adios_global_no_xml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_no_xml_LDADD += $(ADIOSLIB_LDADD)
adios_global_aggregate_by_color_SOURCES = adios_global_aggregate_by_color.c
adios_global_aggregate_by_color_LDADD = $(top_builddir)/src/libadios.a
adios_global_aggregate_by_color_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_aggregate_by_color_LDADD += $(ADIOSLIB_LDADD)
adios_read_global_no_xml_SOURCES = adios_read_global_no_xml.c
adios_read_global_no_xml_LDADD = $(top_builddir)/src/libadiosread.a
adios_read_global_no_xml_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
adios_read_global_no_xml_LDADD += $(ADIOSREADLIB_LDADD)
no_xml_write_byid_SOURCES = no_xml_write_byid.c
no_xml_write_byid_LDADD = $(top_builddir)/src/libadios.a
no_xml_write_byid_LDFLAGS = $(ADIOSLIB_LDFLAGS)
no_xml_write_byid_LDADD += $(ADIOSLIB_LDADD)
read_no_xml_write_byid_SOURCES = read_no_xml_write_byid.c
read_no_xml_write_byid_LDADD = $(top_builddir)/src/libadiosread.a
read_no_xml_write_byid_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
read_no_xml_write_byid_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_global.xml gread_temperature.ch gwrite_temperature.ch job.pbs


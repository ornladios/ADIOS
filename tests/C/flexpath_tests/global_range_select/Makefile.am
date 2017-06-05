AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(srcdir)/../include

AUTOMAKE_OPTIONS = no-dependencies subdir-objects

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/run_tests.pl $(builddir)

noinst_PROGRAMS = writer_range reader_range

writer_range_SOURCES = writer.c ../common/utils.c
writer_range_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
writer_range_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
reader_range_SOURCES = reader.c ../common/utils.c
reader_range_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
reader_range_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = test_config_flex.xml run_tests.pl

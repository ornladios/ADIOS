AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(srcdir)/../include

AUTOMAKE_OPTIONS = no-dependencies subdir-objects

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/readme.txt $(builddir)

noinst_PROGRAMS = writer_two_streams reader_two_streams

writer_two_streams_SOURCES = writer.c ../common/utils.c
writer_two_streams_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
writer_two_streams_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
reader_two_streams_SOURCES = reader.c ../common/utils.c
reader_two_streams_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
reader_two_streams_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)

CLEANFILES = *.bp
CC=$(MPICC)

EXTRA_DIST = cfg.h readme.txt


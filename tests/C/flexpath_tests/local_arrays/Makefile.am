AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(srcdir)/../include

AUTOMAKE_OPTIONS = no-dependencies subdir-objects

noinst_PROGRAMS = writer_local reader_local

writer_local_SOURCES = writer.c ../common/utils.c
writer_local_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
writer_local_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
reader_local_SOURCES = reader.c ../common/utils.c
reader_local_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
reader_local_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)

CLEANFILES = *.bp
CC=$(MPICC)

EXTRA_DIST=cfg.h

AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/job.pbs $(builddir)

noinst_PROGRAMS = stat_write stat_read

stat_write_SOURCES = stat_write.c
stat_write_LDADD = $(top_builddir)/src/libadios.a 
stat_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
stat_write_LDADD += $(ADIOSLIB_LDADD)
stat_read_SOURCES = stat_read.c
stat_read_LDADD = $(top_builddir)/src/libadiosread.a 
stat_read_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
stat_read_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = stat.xml gwrite_stat.ch gread_stat.ch job.pbs

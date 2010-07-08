INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = stat_write stat_read

stat_write_SOURCES = stat_write.c
stat_write_LDADD = $(top_builddir)/src/libadios.a 
stat_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
stat_write_LDADD += $(ADIOSLIB_LDADD)
stat_read_SOURCES = stat_read.c
stat_read_LDADD = $(top_builddir)/src/libadios.a 
stat_read_LDFLAGS = $(ADIOSLIB_LDFLAGS)
stat_read_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = stat.xml gwrite_stat.ch gread_stat.ch job.pbs

AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)
	
noinst_PROGRAMS = scalars_write scalars_read scalars_write_step

scalars_write_SOURCES = scalars_write.c
scalars_write_LDADD = $(top_builddir)/src/libadios.a 
scalars_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
scalars_write_LDADD += $(ADIOSLIB_LDADD)
scalars_read_SOURCES = scalars_read.c
scalars_read_LDADD = $(top_builddir)/src/libadiosread.a 
scalars_read_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
scalars_read_LDADD += $(ADIOSREADLIB_LDADD)
scalars_write_step_SOURCES = scalars_write.c
scalars_write_step_LDADD = $(top_builddir)/src/libadios.a
scalars_write_step_LDFLAGS = $(ADIOSLIB_LDFLAGS)
scalars_write_step_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = scalars.xml gwrite_scalars.ch gread_scalars.ch job.pbs

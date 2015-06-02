AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = scalars_write scalars_read

scalars_write_SOURCES = scalars_write.F90
scalars_write_LDADD = $(top_builddir)/src/libadiosf.a 
scalars_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
scalars_write_LDADD += $(ADIOSLIB_LDADD)

scalars_read_SOURCES = scalars_read.F90
scalars_read_LDADD = $(top_builddir)/src/libadiosreadf.a 
scalars_read_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
scalars_read_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = scalars.xml gwrite_scalars.fh gread_scalars.fh job.pbs

FC=$(MPIFC)

FCLINK = \
        $(LIBTOOL) --mode=link --tag F77 $(FC) \
        $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

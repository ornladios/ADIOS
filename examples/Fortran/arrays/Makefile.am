AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public  -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = arrays_write arrays_read

arrays_write_SOURCES = arrays_write.F90
arrays_write_LDADD = $(top_builddir)/src/libadiosf.a 
arrays_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
arrays_write_LDADD += $(ADIOSLIB_LDADD)

arrays_read_SOURCES = arrays_read.F90
arrays_read_LDADD = $(top_builddir)/src/libadiosreadf.a 
arrays_read_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
arrays_read_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = arrays.xml gwrite_arrays.fh gread_arrays.fh job.pbs

FC=$(MPIFC)

FCLINK = \
        $(LIBTOOL) --mode=link --tag F77 $(FC) \
        $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

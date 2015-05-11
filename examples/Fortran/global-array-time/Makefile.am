AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = adios_globaltime

adios_globaltime_SOURCES = adios_globaltime.F90
adios_globaltime_LDADD = $(top_builddir)/src/libadiosf.a 
adios_globaltime_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_globaltime_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

EXTRA_DIST = adios_globaltime.xml gread_restart.fh gwrite_restart.fh job.pbs

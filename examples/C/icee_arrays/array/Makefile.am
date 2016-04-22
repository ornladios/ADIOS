AM_CPPFLAGS = $(all_includes) -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/README.txt $(builddir)

noinst_PROGRAMS = adios_write adios_read

adios_write_SOURCES = adios_write.c
adios_write_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
adios_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_SOURCES = adios_read.c
adios_read_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
adios_read_LDFLAGS = $(ADIOSLIB_LDFLAGS)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = README.txt

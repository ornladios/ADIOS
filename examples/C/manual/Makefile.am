AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = 1_nonadios_example 2_adios_write 3_adios_read 4_adios_nfiles

1_nonadios_example_SOURCES = 1_nonadios_example.c
1_nonadios_example_LDADD = $(top_builddir)/src/libadios.a 
1_nonadios_example_LDFLAGS = $(ADIOSLIB_LDFLAGS)
1_nonadios_example_LDADD += $(ADIOSLIB_LDADD)
2_adios_write_SOURCES = 2_adios_write.c
2_adios_write_LDADD = $(top_builddir)/src/libadios.a 
2_adios_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
2_adios_write_LDADD += $(ADIOSLIB_LDADD)
3_adios_read_SOURCES = 3_adios_read.c
3_adios_read_LDADD = $(top_builddir)/src/libadios.a 
3_adios_read_LDFLAGS = $(ADIOSLIB_LDFLAGS)
3_adios_read_LDADD += $(ADIOSLIB_LDADD)
4_adios_nfiles_SOURCES = 4_adios_nfiles.c
4_adios_nfiles_LDADD = $(top_builddir)/src/libadios.a 
4_adios_nfiles_LDFLAGS = $(ADIOSLIB_LDFLAGS)
4_adios_nfiles_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = config.xml gread_temperature.ch gwrite_temperature.ch

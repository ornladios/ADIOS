AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = attributes_write attributes_read

attributes_write_SOURCES = attributes_write.c
attributes_write_LDADD = $(top_builddir)/src/libadios.a 
attributes_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
attributes_write_LDADD += $(ADIOSLIB_LDADD)
attributes_read_SOURCES = attributes_read.c
attributes_read_LDADD = $(top_builddir)/src/libadiosread.a
attributes_read_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
attributes_read_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = attributes.xml gwrite_temperature.ch gread_temperature.ch


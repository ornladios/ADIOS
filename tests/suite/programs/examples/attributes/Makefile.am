AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp -p $(srcdir)/*.xml $(builddir)

check_PROGRAMS = attributes_write_C attributes_read_C

attributes_write_C_SOURCES = attributes_write_C.c
attributes_write_C_LDADD = $(top_builddir)/src/libadios.a 
attributes_write_C_LDFLAGS = $(ADIOSLIB_LDFLAGS)
attributes_write_C_LDADD += $(ADIOSLIB_LDADD)
attributes_write_C.o: gwrite_temperature.ch 
gwrite_temperature.ch: attributes_C.xml
	$(top_builddir)/utils/gpp/gpp.py $(srcdir)/attributes_C.xml

attributes_read_C_SOURCES = attributes_read_C.c
attributes_read_C_LDADD = $(top_builddir)/src/libadiosread.a
attributes_read_C_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
attributes_read_C_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp *.ch *.fh
CC=$(MPICC)
EXTRA_DIST = attributes_C.xml 


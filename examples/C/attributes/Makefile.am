INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = attributes_write attributes_read

attributes_write_SOURCES = attributes_write.c
attributes_write_LDADD = $(top_builddir)/src/libadios.a 
attributes_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
attributes_write_LDADD += $(ADIOSLIB_LDADD)
attributes_read_SOURCES = attributes_read.c
attributes_read_LDADD = $(top_builddir)/src/libadiosread.a
attributes_read_LDFLAGS = $(ADIOSLIB_LDFLAGS)
attributes_read_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = attributes.xml gwrite_temperature.ch gread_temperature.ch


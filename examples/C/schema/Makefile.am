INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I$(top_builddir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = schema_write unstructured

schema_write_SOURCES = schema_write.c
schema_write_LDADD = $(top_builddir)/src/libadios.a 
schema_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
schema_write_LDADD += $(ADIOSLIB_LDADD)

unstructured_SOURCES = unstructured.c
unstructured_LDADD = $(top_builddir)/src/libadios.a 
unstructured_LDFLAGS = $(ADIOSLIB_LDFLAGS)
unstructured_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = uniform.xml structured.xml rectilinear.xml unstructured.xml gwrite_temperature.ch


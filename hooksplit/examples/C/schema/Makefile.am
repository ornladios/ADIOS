INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I$(top_builddir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = uniform rectilinear structured unstructured

uniform_SOURCES = uniform.c
uniform_LDADD = $(top_builddir)/src/libadios.a 
uniform_LDFLAGS = $(ADIOSLIB_LDFLAGS)
uniform_LDADD += $(ADIOSLIB_LDADD)

rectilinear_SOURCES = rectilinear.c
rectilinear_LDADD = $(top_builddir)/src/libadios.a 
rectilinear_LDFLAGS = $(ADIOSLIB_LDFLAGS)
rectilinear_LDADD += $(ADIOSLIB_LDADD)

structured_SOURCES = structured.c
structured_LDADD = $(top_builddir)/src/libadios.a 
structured_LDFLAGS = $(ADIOSLIB_LDFLAGS)
structured_LDADD += $(ADIOSLIB_LDADD)

unstructured_SOURCES = unstructured.c
unstructured_LDADD = $(top_builddir)/src/libadios.a 
unstructured_LDFLAGS = $(ADIOSLIB_LDFLAGS)
unstructured_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = uniform.xml structured.xml rectilinear.xml unstructured.xml 


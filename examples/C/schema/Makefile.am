INCLUDES = $(all_includes)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = uniform rectilinear structured unstructured triangle2d

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

triangle2d_SOURCES = triangle2d.c
triangle2d_LDADD = $(top_builddir)/src/libadios.a 
triangle2d_LDFLAGS = $(ADIOSLIB_LDFLAGS)
triangle2d_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = uniform.xml structured.xml rectilinear.xml unstructured.xml \
             triangle2d.xml triangle2d.readme triangle2d.n4.png


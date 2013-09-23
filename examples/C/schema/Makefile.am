INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I$(top_builddir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = uniform rectilinear structured unstructured uniform2d rectilinear2d structured2d tri2d triangle2d

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

uniform2d_SOURCES = uniform2d.c
uniform2d_LDADD = $(top_builddir)/src/libadios.a
uniform2d_LDFLAGS = $(ADIOSLIB_LDFLAGS)
uniform2d_LDADD += $(ADIOSLIB_LDADD)

rectilinear2d_SOURCES =  rectilinear2d.c
rectilinear2d_LDADD = $(top_builddir)/src/libadios.a
rectilinear2d_LDFLAGS = $(ADIOSLIB_LDFLAGS)
rectilinear2d_LDADD += $(ADIOSLIB_LDADD)

structured2d_SOURCES =  structured2d.c
structured2d_LDADD = $(top_builddir)/src/libadios.a
structured2d_LDFLAGS = $(ADIOSLIB_LDFLAGS)
structured2d_LDADD += $(ADIOSLIB_LDADD)

tri2d_SOURCES = tri2d.c
tri2d_LDADD = $(top_builddir)/src/libadios.a
tri2d_LDFLAGS = $(ADIOSLIB_LDFLAGS)
tri2d_LDADD += $(ADIOSLIB_LDADD)

triangle2d_SOURCES = triangle2d.c
triangle2d_LDADD = $(top_builddir)/src/libadios.a 
triangle2d_LDFLAGS = $(ADIOSLIB_LDFLAGS)
triangle2d_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = uniform.xml structured.xml rectilinear.xml unstructured.xml \
			 unifrom2d.xml unifrom2d.readme unifrom2d.n4.png \
			 rectilinear2d.xml rectilinear2d.readme rectilinear2d.png \
			 structured2d.xml structured2d.readme \
			 tri2d.xml tri2d.readme tri2d.n.png \
             triangle2d.xml triangle2d.readme triangle2d.n4.png


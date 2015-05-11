AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = uniform2d rectilinear2d structured2d tri2d triangle2d uniform2d_noxml rectilinear2d_noxml structured2d_noxml tri2d_noxml

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

uniform2d_noxml_LDADD = $(top_builddir)/src/libadios.a
uniform2d_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
uniform2d_noxml_LDADD += $(ADIOSLIB_LDADD)

rectilinear2d_noxml_LDADD = $(top_builddir)/src/libadios.a
rectilinear2d_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
rectilinear2d_noxml_LDADD += $(ADIOSLIB_LDADD)

structured2d_noxml_LDADD = $(top_builddir)/src/libadios.a
structured2d_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
structured2d_noxml_LDADD += $(ADIOSLIB_LDADD)

tri2d_noxml_LDADD = $(top_builddir)/src/libadios.a
tri2d_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
tri2d_noxml_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = uniform2d.xml uniform2d.readme uniform2d.png \
	rectilinear2d.xml rectilinear2d.readme rectilinear2d.png \
	structured2d.xml structured2d.readme \
	tri2d.xml tri2d.readme tri2d.n.png \
	triangle2d.xml triangle2d.readme triangle2d.n4.png


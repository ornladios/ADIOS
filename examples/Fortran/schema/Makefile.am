AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

#all-local:
#	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = uniform2d_f_noxml rectilinear2d_f_noxml structured2d_f_noxml tri2d_f_noxml tri2d_f_noxml_seperate 

uniform2d_f_noxml_SOURCES = uniform2d_noxml.F90
uniform2d_f_noxml_LDADD = $(top_builddir)/src/libadiosf.a
uniform2d_f_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
uniform2d_f_noxml_LDADD += $(ADIOSLIB_LDADD)

rectilinear2d_f_noxml_SOURCES = rectilinear2d_noxml.F90
rectilinear2d_f_noxml_LDADD = $(top_builddir)/src/libadiosf.a
rectilinear2d_f_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
rectilinear2d_f_noxml_LDADD += $(ADIOSLIB_LDADD)

structured2d_f_noxml_SOURCES = structured2d_noxml.F90
structured2d_f_noxml_LDADD = $(top_builddir)/src/libadiosf.a
structured2d_f_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
structured2d_f_noxml_LDADD += $(ADIOSLIB_LDADD)

tri2d_f_noxml_SOURCES = tri2d_noxml.F90 
tri2d_f_noxml_LDADD = $(top_builddir)/src/libadiosf.a
tri2d_f_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
tri2d_f_noxml_LDADD += $(ADIOSLIB_LDADD)

tri2d_f_noxml_seperate_SOURCES = tri2d_noxml_seperate.F90
tri2d_f_noxml_seperate_LDADD = $(top_builddir)/src/libadiosf.a
tri2d_f_noxml_seperate_LDFLAGS = $(ADIOSLIB_LDFLAGS)
tri2d_f_noxml_seperate_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
FC=$(MPIFC)

FCLINK = \
    $(LIBTOOL) --mode=link --tag F77 $(FC) \
    $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@


INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

#all-local:
#	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS = tri2d_f_noxml 

tri2d_f_noxml_SOURCES = tri2d_noxml.F90 
tri2d_f_noxml_LDADD = $(top_builddir)/src/libadiosf.a
tri2d_f_noxml_LDFLAGS = $(ADIOSLIB_LDFLAGS)
tri2d_f_noxml_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
FC=$(MPIFC)

FCLINK = \
    $(LIBTOOL) --mode=link --tag F77 $(FC) \
    $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@


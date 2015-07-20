AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp -p $(srcdir)/*.xml $(builddir)
	
check_PROGRAMS = scalars_write_C scalars_read_C

scalars_write_C_SOURCES = scalars_write_C.c 
scalars_write_C_LDADD = $(top_builddir)/src/libadios.a 
scalars_write_C_LDFLAGS = $(ADIOSLIB_LDFLAGS)
scalars_write_C_LDADD += $(ADIOSLIB_LDADD)
scalars_write_C.o: gwrite_scalars.ch 
gwrite_scalars.ch: scalars_C.xml
	$(top_builddir)/utils/gpp/gpp.py $(srcdir)/scalars_C.xml

scalars_read_C_SOURCES = scalars_read_C.c
scalars_read_C_LDADD = $(top_builddir)/src/libadiosread.a 
scalars_read_C_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
scalars_read_C_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp *.ch *.fh
CC=$(MPICC)
EXTRA_DIST = scalars_C.xml scalars_F.xml 


if BUILD_FORTRAN
check_PROGRAMS += scalars_write_F scalars_read_F
scalars_write_F_SOURCES = scalars_write_F.F90
scalars_write_F_LDADD = $(top_builddir)/src/libadiosf.a 
scalars_write_F_LDFLAGS = $(ADIOSLIB_LDFLAGS)
scalars_write_F_LDADD += $(ADIOSLIB_LDADD)
scalars_write_F.o: gwrite_scalars.fh 
gwrite_scalars.fh: scalars_F.xml
	$(top_builddir)/utils/gpp/gpp.py $(srcdir)/scalars_F.xml

scalars_read_F_SOURCES = scalars_read_F.F90
scalars_read_F_LDADD = $(top_builddir)/src/libadiosreadf.a 
scalars_read_F_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
scalars_read_F_LDADD += $(ADIOSREADLIB_LDADD)

endif

FC=$(MPIFC)
FCLINK = \
        $(LIBTOOL) --mode=link --tag F77 $(FC) \
        $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@




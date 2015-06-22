AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public  -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp -p $(srcdir)/*.xml $(builddir)

check_PROGRAMS = global_array_write_C global_array_read_C 

global_array_write_C_SOURCES = global_array_write_C.c
global_array_write_C_LDADD = $(top_builddir)/src/libadios.a 
global_array_write_C_LDFLAGS = $(ADIOSLIB_LDFLAGS)
global_array_write_C_LDADD += $(ADIOSLIB_LDADD)
global_array_write_C.o: gwrite_temperature.ch 
gwrite_temperature.ch: global_array_C.xml
	$(top_builddir)/utils/gpp/gpp.py $(srcdir)/global_array_C.xml

global_array_read_C_SOURCES = global_array_read_C.c
global_array_read_C_LDADD = $(top_builddir)/src/libadiosread.a 
global_array_read_C_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
global_array_read_C_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp *.ch *.fh
CC=$(MPICC)
EXTRA_DIST = global_array_C.xml global_array_F.xml 


if BUILD_FORTRAN
check_PROGRAMS += global_array_write_F 

global_array_write_F_SOURCES = global_array_write_F.F90
global_array_write_F_LDADD = $(top_builddir)/src/libadiosf.a 
global_array_write_F_LDFLAGS = $(ADIOSLIB_LDFLAGS)
global_array_write_F_LDADD += $(ADIOSLIB_LDADD)
global_array_write_F.o: gwrite_temperature.fh 
gwrite_temperature.fh: global_array_F.xml
	$(top_builddir)/utils/gpp/gpp.py $(srcdir)/global_array_F.xml

endif

FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@



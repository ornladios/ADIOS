INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src #-I/ccs/home/sklasky/parallel-netcdf/jaguar/include

AUTOMAKE_OPTIONS = no-dependencies

if BUILD_FORTRAN
noinst_PROGRAMS=genbp bp_read_f bp_read_c
else
noinst_PROGRAMS=genbp bp_read_c
endif

# read_coords read_phi read_parts_global read_gena read_gena_pncd read_pixie_all

genbp_SOURCES = genbp.c
genbp_CPPFLAGS = $(ADIOSLIB_CPPFLAGS)
genbp_CFLAGS = $(ADIOSLIB_CFLAGS)
genbp_LDADD = $(top_builddir)/src/libadios.a $(top_builddir)/src/libadiosread.a $(ADIOSLIB_LDADD)
genbp_LDFLAGS = $(ADIOSLIB_LDFLAGS) 

if BUILD_FORTRAN
bp_read_f_SOURCES = bp_read_f.F90
bp_read_f_CPPFLAGS = $(ADIOSLIB_CPPFLAGS)
bp_read_f_CFLAGS = $(ADIOSLIB_CFLAGS)
bp_read_f_LDADD = $(top_builddir)/src/libadiosreadf.a $(ADIOSREADLIB_LDADD)
bp_read_f_LDFLAGS = $(ADIOSREADLIB_LDFLAGS) 
endif

bp_read_c_SOURCES = bp_read_c.c
bp_read_c_CPPFLAGS = $(ADIOSLIB_CPPFLAGS)
bp_read_c_CFLAGS = $(ADIOSLIB_CFLAGS)
bp_read_c_LDADD = $(top_builddir)/src/libadiosread.a $(ADIOSREADLIB_LDADD)
bp_read_c_LDFLAGS = $(ADIOSREADLIB_LDFLAGS) 

FC=$(MPIFC)
CC=$(MPICC)

FCLINK = \
        $(LIBTOOL) --mode=link --tag F77 $(FC) \
        $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

if BUILD_FORTRAN
CLEANFILES = *.o bp_read_c bp_read_f genbp read_coords read_phi read_parts_global
else
CLEANFILES = *.o bp_read_c genbp read_coords read_phi read_parts_global
endif


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
genbp_LDFLAGS = $(MXML_LDFLAGS)
genbp_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_NETCDF
genbp_LDADD +=  $(NETCDF_LIBS)
genbp_LDFLAGS += $(NETCDF_LDFLAGS)
endif
if HAVE_NSSI
genbp_LDADD +=  $(NSSI_LIBS) $(EXTRA_CXX_LIBS)
genbp_LDFLAGS += $(NSSI_LDFLAGS)
endif
if HAVE_PHDF5
genbp_LDADD +=  $(PHDF5_LIBS)
genbp_LDFLAGS += $(PHDF5_LDFLAGS)
endif
if HAVE_DART
genbp_CPPFLAGS = $(DART_CPPFLAGS)
genbp_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
genbp_LDADD += $(DART_LIBS)
endif

if BUILD_FORTRAN
bp_read_f_SOURCES = bp_read_f.F90
bp_read_f_LDFLAGS = 
bp_read_f_LDADD = $(top_builddir)/src/libadiosreadf.a
if HAVE_DMALLOC
bp_read_f_CFLAGS = -DDMALLOC  $(DMALLOC_CPPFLAGS)
bp_read_f_LDFLAGS +=  $(DMALLOC_LDFLAGS)
bp_read_f_LDADD += -ldmalloc
endif
if HAVE_DART
bp_read_f_CPPFLAGS = $(DART_CPPFLAGS)
bp_read_f_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
bp_read_f_LDADD += $(DART_LIBS)
endif
endif

#read_coords_SOURCES = read_coords.F90
#read_coords_LDADD =  $(top_builddir)/src/libadiosread.a
#
#read_phi_SOURCES = read_phi.F90
#read_phi_LDADD = $(top_builddir)/src/libadiosread.a
#
#read_pixie_all_SOURCES = read_pixie_all.F90
#read_pixie_all_LDADD = $(top_builddir)/src/libadiosread.a

#read_parts_global_SOURCES = read_parts_global.F90
#read_parts_global_LDFLAGS = $(MXML_LDFLAGS)
#read_parts_global_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
#if HAVE_PHDF5
#read_parts_global_LDADD +=  $(PHDF5_LIBS)
#read_parts_global_LDFLAGS += $(PHDF5_LDFLAGS)
#endif
#
#read_gena_SOURCES = read_gena.F90
#read_gena_LDFLAGS = $(MXML_LDFLAGS)
#read_gena_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
#if HAVE_PHDF5
#read_gena_LDADD +=  $(PHDF5_LIBS)
#read_gena_LDFLAGS += $(PHDF5_LDFLAGS)
#endif

#read_gena_pncd_SOURCES = read_gena_pncd.c
#read_gena_pncd_LDFLAGS = $(MXML_LDFLAGS)
#read_gena_pncd_LDADD = /ccs/home/sklasky/parallel-netcdf/jaguar/lib/libpnetcdf.a 
#if HAVE_PHDF5
#read_gena_pncd_LDADD +=  $(PHDF5_LIBS)
#read_gena_pncd_LDFLAGS += $(PHDF5_LDFLAGS)
#endif


bp_read_c_SOURCES = bp_read_c.c
bp_read_c_LDFLAGS = 
bp_read_c_LDADD = $(top_builddir)/src/libadiosread.a
if HAVE_DMALLOC
bp_read_c_CFLAGS = -DDMALLOC  $(DMALLOC_CPPFLAGS)
bp_read_c_LDFLAGS +=  $(DMALLOC_LDFLAGS)
bp_read_c_LDADD += -ldmalloc
endif
if HAVE_DART
bp_read_c_CPPFLAGS = $(DART_CPPFLAGS)
bp_read_c_LDFLAGS += $(DART_CPPFLAGS) $(DART_LDFLAGS)
bp_read_c_LDADD += $(DART_LIBS)
endif

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


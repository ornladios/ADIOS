INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I/ccs/home/sklasky/parallel-netcdf/jaguar/include

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=genbp bp_read_f bp_read_c read_coords read_phi read_parts_global read_gena read_gena_pncd read_pixie_all

genbp_SOURCES = genbp.c
genbp_LDFLAGS = $(MXML_LDFLAGS)
genbp_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
genbp_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
genbp_LDFLAGS += $(PHDF5_LDFLAGS)
endif

bp_read_f_SOURCES = bp_read_f.F90
bp_read_f_LDFLAGS = $(MXML_LDFLAGS)
bp_read_f_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
bp_read_f_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
bp_read_f_LDFLAGS += $(PHDF5_LDFLAGS)
endif

read_coords_SOURCES = read_coords.F90
read_coords_LDFLAGS = $(MXML_LDFLAGS)
read_coords_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
read_coords_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
read_coords_LDFLAGS += $(PHDF5_LDFLAGS)
endif

read_phi_SOURCES = read_phi.F90
read_phi_LDFLAGS = $(MXML_LDFLAGS)
read_phi_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
read_phi_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
read_phi_LDFLAGS += $(PHDF5_LDFLAGS)
endif

read_pixie_all_SOURCES = read_pixie_all.F90
read_pixie_all_LDFLAGS = $(MXML_LDFLAGS)
read_pixie_all_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
read_pixie_all_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
read_pixie_all_LDFLAGS += $(PHDF5_LDFLAGS)
endif

read_parts_global_SOURCES = read_parts_global.F90
read_parts_global_LDFLAGS = $(MXML_LDFLAGS)
read_parts_global_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
read_parts_global_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
read_parts_global_LDFLAGS += $(PHDF5_LDFLAGS)
endif

read_gena_SOURCES = read_gena.F90
read_gena_LDFLAGS = $(MXML_LDFLAGS)
read_gena_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
read_gena_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
read_gena_LDFLAGS += $(PHDF5_LDFLAGS)
endif

read_gena_pncd_SOURCES = read_gena_pncd.c
read_gena_pncd_LDFLAGS = $(MXML_LDFLAGS)
read_gena_pncd_LDADD = /ccs/home/sklasky/parallel-netcdf/jaguar/lib/libpnetcdf.a 
if HAVE_PHDF5
read_gena_pncd_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
read_gena_pncd_LDFLAGS += $(PHDF5_LDFLAGS)
endif


bp_read_c_SOURCES = bp_read_c.c
bp_read_c_LDFLAGS = $(MXML_LDFLAGS)
bp_read_c_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
bp_read_c_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
bp_read_c_LDFLAGS += $(PHDF5_LDFLAGS)
endif

FC=$(MPIFC)
CC=$(MPICC)

FCLINK = \
        $(LIBTOOL) --mode=link --tag F77 $(FC) \
        $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

CLEANFILES = *.o bp_read_c bp_read_f genbp read_coords read_phi read_parts_global

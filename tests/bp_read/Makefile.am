INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=genbp bp_read_f bp_read_c

genbp_SOURCES = genbp.c
genbp_LDFLAGS = $(MXML_LDFLAGS)
genbp_LDADD = $(top_builddir)/src/libadios.a $(MXML_LIBS) $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
genbp_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
genbp_LDFLAGS += $(PHDF5_LDFLAGS)
endif

bp_read_f_SOURCES = bp_read_f.F90
bp_read_f_LDFLAGS = $(MXML_LDFLAGS)
bp_read_f_LDADD = $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
bp_read_f_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
bp_read_f_LDFLAGS += $(PHDF5_LDFLAGS)
endif

bp_read_c_SOURCES = bp_read_c.c
bp_read_c_LDFLAGS = $(MXML_LDFLAGS)
bp_read_c_LDADD = $(top_builddir)/src/libadiosread.a
if HAVE_PHDF5
bp_read_c_LDADD += $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
bp_read_c_LDFLAGS += $(PHDF5_LDFLAGS)
endif

FC=$(MPIFC)

FCLINK = \
        $(LIBTOOL) --mode=link --tag F77 $(FC) \
        $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

CLEANFILES = *.o bp_read_c bp_read_f genbp

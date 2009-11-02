INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = 1_nonadios_example 2_adios_posix 3_adios_nfiles

1_nonadios_example_SOURCES = 1_nonadios_example.c
1_nonadios_example_LDADD = $(top_builddir)/src/libadios.a 
1_nonadios_example_LDFLAGS = $(MXML_LDFLAGS)
1_nonadios_example_LDADD += $(MXML_LIBS)
2_adios_posix_SOURCES = 2_adios_posix.c
2_adios_posix_LDADD = $(top_builddir)/src/libadios.a 
2_adios_posix_LDFLAGS = $(MXML_LDFLAGS)
2_adios_posix_LDADD += $(MXML_LIBS)
3_adios_nfiles_SOURCES = 3_adios_nfiles.c
3_adios_nfiles_LDADD = $(top_builddir)/src/libadios.a 
3_adios_nfiles_LDFLAGS = $(MXML_LDFLAGS)
3_adios_nfiles_LDADD += $(MXML_LIBS)
if HAVE_PHDF5
1_nonadios_example_LDFLAGS += $(PHDF5_LDFLAGS)
1_nonadios_example_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
2_adios_posix_LDFLAGS += $(PHDF5_LDFLAGS)
2_adios_posix_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
3_adios_nfiles_LDFLAGS += $(PHDF5_LDFLAGS)
3_adios_nfiles_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)


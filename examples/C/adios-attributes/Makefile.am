INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = adios_attributes

adios_attributes_SOURCES = adios_attributes.c
adios_attributes_LDADD = $(top_builddir)/src/libadios.a 
adios_attributes_LDFLAGS = $(MXML_LDFLAGS)
adios_attributes_LDADD += $(MXML_LIBS)
if HAVE_PHDF5
adios_attributes_LDFLAGS += $(PHDF5_LDFLAGS)
adios_attributes_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)


INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = attributes_write attributes_read

attributes_write_SOURCES = attributes_write.c
attributes_write_LDADD = $(top_builddir)/src/libadios.a 
attributes_write_LDFLAGS = $(MXML_LDFLAGS)
attributes_write_LDADD += $(MXML_LIBS)
attributes_read_SOURCES = attributes_read.c
attributes_read_LDADD = $(top_builddir)/src/libadiosread.a
attributes_read_LDFLAGS = $(MXML_LDFLAGS)
attributes_read_LDADD += $(MXML_LIBS)

if HAVE_PHDF5
attributes_write_LDFLAGS += $(PHDF5_LDFLAGS)
attributes_write_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
attributes_read_LDFLAGS += $(PHDF5_LDFLAGS)
attributes_read_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = attributes.xml 


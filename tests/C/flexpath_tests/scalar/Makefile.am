INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I$(top_builddir)/src/public -I$(srcdir)/../include

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = writer reader

writer_SOURCES = writer.c ../common/utils.c
writer_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
writer_LDFLAGS = $(ADIOSLIB_LDFLAGS)
reader_SOURCES = reader.c ../common/utils.c
reader_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
reader_LDFLAGS = $(ADIOSLIB_LDFLAGS)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = arrays.xml 

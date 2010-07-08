INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = scalars_write scalars_read

scalars_write_SOURCES = scalars_write.c
scalars_write_LDADD = $(top_builddir)/src/libadios.a 
scalars_write_LDFLAGS = $(ADIOSLIB_LDFLAGS)
scalars_write_LDADD += $(ADIOSLIB_LDADD)
scalars_read_SOURCES = scalars_read.c
scalars_read_LDADD = $(top_builddir)/src/libadios.a 
scalars_read_LDFLAGS = $(ADIOSLIB_LDFLAGS)
scalars_read_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = scalars.xml gwrite_scalars.ch gread_scalars.ch job.pbs

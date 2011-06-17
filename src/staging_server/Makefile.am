INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I$(top_builddir)/src/ardma

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=bp_staging_server

bp_staging_server_SOURCES = options.c mem.c precedence.c transport_thread.c bpworker_thread.c bpwriter_thread.c bpindexing_thread.c queue.c dirutil.c main.c

bp_staging_server_CPPFLAGS = $(ADIOSLIB_CPPFLAGS)
bp_staging_server_CFLAGS = $(ADIOSLIB_CFLAGS)
bp_staging_server_LDADD = $(top_builddir)/src/libadios.a $(top_builddir)/src/ardma/libardma_server.a $(ADIOSLIB_LDADD)
bp_staging_server_LDFLAGS = $(ADIOSLIB_LDFLAGS) -L$(top_builddir)/src/ardma

#CLEANFILES = restart.bp
CC=$(MPICC)

#EXTRA_DIST = config_c.xml

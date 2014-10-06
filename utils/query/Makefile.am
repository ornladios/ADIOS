INCLUDES = $(all_includes)
#

INCLUDES += -I$(top_srcdir)/src/public ${ADIOSREADLIB_CPPFLAGS}

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = index_fastbit

index_fastbit_SOURCES = indexByFastbit.c #fastbit_adios.c fastbit_adios.h

index_fastbit_CPPFLAGS = $(ADIOSLIB_CPPFLAGS) $(ADIOSLIB_CFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS)
index_fastbit_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
index_fastbit_LDFLAGS = $(ADIOSLIB_LDFLAGS)




CLEANFILES = *.idx
CC=$(MPICC)

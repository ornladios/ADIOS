INCLUDES = $(all_includes)
#

INCLUDES += -I$(top_srcdir)/src/public ${ADIOSREADLIB_CPPFLAGS}

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = index_fastbit

index_fastbit_SOURCES = indexByFastbit.c #fastbit_adios.c fastbit_adios.h

index_fastbit_LDADD = $(top_builddir)/src/libadiosread.a $(top_builddir)/src/libadios.a  -lmxml -lfastbit
index_fastbit_LDFLAGS = $(ADIOSREADLIB_LDFLAGS) $(ADIOSLIB_LDFLAGS)




CLEANFILES = *.idx
CC=$(MPICC)
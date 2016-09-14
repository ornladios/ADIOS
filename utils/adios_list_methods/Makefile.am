AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public  -I$(top_srcdir)/src -I$(top_srcdir)/src/core

AUTOMAKE_OPTIONS = no-dependencies

#
# Executables for the sequential libraries
#
bin_PROGRAMS = adios_list_methods_readonly_nompi
adios_list_methods_readonly_nompi_SOURCES = adios_list_methods.c
adios_list_methods_readonly_nompi_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSREADLIB_SEQ_CPPFLAGS) $(ADIOSREADLIB_SEQ_CFLAGS) 
adios_list_methods_readonly_nompi_LDFLAGS = $(ADIOSREADLIB_SEQ_LDFLAGS) 
adios_list_methods_readonly_nompi_LDADD =  $(top_builddir)/src/libadiosread_nompi.a $(ADIOSREADLIB_SEQ_LDADD)

if BUILD_WRITE
bin_PROGRAMS += adios_list_methods_nompi 
adios_list_methods_nompi_SOURCES = adios_list_methods.c
adios_list_methods_nompi_CPPFLAGS = $(AM_CPPFLAGS) ${MACRODEFFLAG}WRITE $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) 
adios_list_methods_nompi_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
adios_list_methods_nompi_LDADD =  $(top_builddir)/src/libadios_nompi.a $(ADIOSLIB_SEQ_LDADD)
endif BUILD_WRITE

#
# Executables for the parallel libraries
#
if HAVE_MPI
bin_PROGRAMS += adios_list_methods_readonly
adios_list_methods_readonly_SOURCES = adios_list_methods.c
adios_list_methods_readonly_CPPFLAGS = $(AM_CPPFLAGS) $(ADIOSREADLIB_CPPFLAGS) $(ADIOSREADLIB_CFLAGS)
adios_list_methods_readonly_LDFLAGS = $(ADIOSREADLIB_LDFLAGS) 
adios_list_methods_readonly_LDADD =  $(top_builddir)/src/libadiosread.a $(ADIOSREADLIB_LDADD)

if BUILD_WRITE
bin_PROGRAMS += adios_list_methods adios_list_methods_nompi 
adios_list_methods_SOURCES = adios_list_methods.c
adios_list_methods_CPPFLAGS = $(AM_CPPFLAGS) ${MACRODEFFLAG}WRITE $(ADIOSLIB_CPPFLAGS) $(ADIOSLIB_CFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) 
adios_list_methods_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
adios_list_methods_LDADD =  $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
endif BUILD_WRITE
override CC=$(MPICC)
override CXX=$(MPICXX)
endif HAVE_MPI



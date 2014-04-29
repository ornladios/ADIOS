INCLUDES = $(all_includes)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/core

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = list_methods list_methods_nompi list_methods_readonly list_methods_readonly_nompi

if BUILD_WRITE
list_methods_SOURCES = list_methods.c
list_methods_CPPFLAGS = ${MACRODEFFLAG}WRITE $(ADIOSLIB_CPPFLAGS) $(ADIOSLIB_CFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) 
list_methods_LDFLAGS = $(ADIOSLIB_LDFLAGS) 
list_methods_LDADD =  $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)

list_methods_nompi_SOURCES = list_methods.c
list_methods_nompi_CPPFLAGS = ${MACRODEFFLAG}WRITE $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS) $(ADIOSLIB_EXTRA_CPPFLAGS) 
list_methods_nompi_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS) 
list_methods_nompi_LDADD =  $(top_builddir)/src/libadios_nompi.a $(ADIOSLIB_SEQ_LDADD)
endif

list_methods_readonly_SOURCES = list_methods.c
list_methods_readonly_CPPFLAGS = $(ADIOSREADLIB_CPPFLAGS) $(ADIOSREADLIB_CFLAGS)
list_methods_readonly_LDFLAGS = $(ADIOSREADLIB_LDFLAGS) 
list_methods_readonly_LDADD =  $(top_builddir)/src/libadiosread.a $(ADIOSREADLIB_LDADD)

list_methods_readonly_nompi_SOURCES = list_methods.c
list_methods_readonly_nompi_CPPFLAGS = $(ADIOSREADLIB_SEQ_CPPFLAGS) $(ADIOSREADLIB_SEQ_CFLAGS) 
list_methods_readonly_nompi_LDFLAGS = $(ADIOSREADLIB_SEQ_LDFLAGS) 
list_methods_readonly_nompi_LDADD =  $(top_builddir)/src/libadiosread_nompi.a $(ADIOSREADLIB_SEQ_LDADD)

#if USE_PARALLEL_COMPILER
CC=$(MPICC)
CXX=$(MPICXX)
#endif



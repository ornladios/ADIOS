INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src -I$(top_builddir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = read_all

read_all_SOURCES = read_all.c
read_all_CPPFLAGS = ${MACRODEFFLAG}_NOMPI ${MACRODEFFLAG}ADIOS_USE_READ_API_1
read_all_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
read_all_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
read_all_LDADD += $(ADIOSREADLIB_LDADD)

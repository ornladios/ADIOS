AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = read_all read_all_v1

read_all_SOURCES = read_all.c
read_all_CPPFLAGS = ${MACRODEFFLAG}_NOMPI $(AM_CPPFLAGS)
read_all_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
read_all_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
read_all_LDADD += $(ADIOSREADLIB_LDADD)

read_all_v1_SOURCES = read_all_v1.c
read_all_v1_CPPFLAGS = ${MACRODEFFLAG}_NOMPI ${MACRODEFFLAG}ADIOS_USE_READ_API_1 $(AM_CPPFLAGS)
read_all_v1_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
read_all_v1_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
read_all_v1_LDADD += $(ADIOSREADLIB_LDADD)

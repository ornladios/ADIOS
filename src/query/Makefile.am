INCLUDES = $(all_includes)
#
# ${ADIOSREADLIB_CPPFLAGS} has all the flags like -DFASTBIT -DALACRITY
# as well as related includes, so it is needed here
# 
INCLUDES += -I$(top_srcdir)/src/public -I$(top_srcdir)/src/ -I$(top_srcdir)/src/${ADIOSREADLIB_CPPFLAGS}

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adiosQuery

adiosQuery_SOURCES = adios_query.c  adios_query_hooks.c  common_query.c  fastbit_adios.c  mock_alacrity.c  myTestCases.c  query_fastbit.c
query_fastbit_SOURCES = fastbit_adios.c query_fastbit.c common_query.c query.c adios_query_hooks.c myTestCases.c
index_fastbit_SOURCES = useindex.c fastbit_adios.c
mock_alacrity_SOURCES = fastbit_adios.c mock_alacrity.c common_query.c query.c adios_query_hooks.c myTestCases.c
adiosQuery_CPPFLAGS = ${MACRODEFFLAG}_NOMPI 
adiosQuery_LDADD = $(top_builddir)/src/libadiosread_nompi.a $(ADIOSREADLIB_LDADD)
adiosQuery_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)

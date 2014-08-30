INCLUDES = $(all_includes)
#
# ${ADIOSREADLIB_CPPFLAGS} has all the flags like -DFASTBIT -DALACRITY
# as well as related includes, so it is needed here
# 
INCLUDES += -I$(top_srcdir)/src/

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = 

# Only compile these test programs if we have FastBit, since they rely on FastBit
if HAVE_FASTBIT
noinst_PROGRAMS += query_fastbit

query_fastbit_SOURCES = ../testCases.c
query_fastbit_CPPFLAGS = ${MACRODEFFLAG}_NOMPI 
query_fastbit_LDADD = $(top_builddir)/src/libadiosread_nompi.a $(ADIOSREADLIB_LDADD)
query_fastbit_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)

endif # HAVE_FASTBIT
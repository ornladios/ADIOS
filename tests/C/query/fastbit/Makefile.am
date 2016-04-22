AM_CPPFLAGS = $(all_includes)
#
# ${ADIOSREADLIB_CPPFLAGS} has all the flags like -DFASTBIT -DALACRITY
# as well as related includes, so it is needed here
# 
AM_CPPFLAGS += -I$(top_builddir)/src/public  -I$(top_srcdir)/src/ -I$(top_srcdir)/src/public/ $(MXML_CPPFLAGS)

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = 

# Only compile these test programs if we have FastBit, since they rely on FastBit
if HAVE_FASTBIT

noinst_PROGRAMS += fastbit_tests

fastbit_tests_SOURCES = fastbit_tests.c
fastbit_tests_CPPFLAGS = ${MACRODEFFLAG}_NOMPI $(AM_CPPFLAGS)
fastbit_tests_LDADD = $(top_builddir)/src/libadiosread_nompi.a $(ADIOSREADLIB_SEQ_LDADD) ${MXML_LIBS}
fastbit_tests_LDFLAGS = $(ADIOSREADLIB_SEQ_LDFLAGS) ${MXML_LDFLAGS}

endif # HAVE_FASTBIT

AM_CPPFLAGS = $(all_includes)
#
# ${ADIOSREADLIB_CPPFLAGS} has all the flags like -DFASTBIT -DALACRITY
# as well as related includes, so it is needed here
# 
AM_CPPFLAGS += -I$(top_builddir)/src/public  -I$(top_srcdir)/src/ -I$(top_srcdir)/src/public/ 

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = 

# Only compile these test programs if we have FastBit, since they rely on FastBit
if HAVE_FASTBIT

noinst_PROGRAMS += fastbit_tests

fastbit_tests_SOURCES = fastbit_tests.c
fastbit_tests_CPPFLAGS = $(AM_CPPFLAGS) ${ADIOSLIB_SEQ_CPPFLAGS} $(ADIOSLIB_EXTRA_CPPFLAGS)
fastbit_tests_LDADD = $(top_builddir)/src/libadios_nompi.a $(ADIOSLIB_SEQ_LDADD) 
fastbit_tests_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)

EXTRA_DIST = query.xml

endif # HAVE_FASTBIT

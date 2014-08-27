INCLUDES = $(all_includes)
INCLUDES += $(ADIOSLIB_CPPFLAGS)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/public -I$(srcdir)/../include

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/readme.txt $(builddir)

noinst_PROGRAMS =

if HAVE_ALACRITY
noinst_PROGRAMS += query_test
query_test_SOURCES = test_query_alac_internal.c 
query_test_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD) 
query_test_LDFLAGS = $(ADIOSLIB_LDFLAGS)  
endif # HAVE_ALACRITY

CLEANFILES = *.bp
CC=$(MPICC)
#EXTRA_DIST = test_config_flex.xml 

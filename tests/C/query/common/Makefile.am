INCLUDES = $(all_includes)
INCLUDES += $(ADIOSLIB_CPPFLAGS)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = build_indexed_dataset

build_indexed_dataset_SOURCES = build_indexed_dataset.c 
build_indexed_dataset_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD) 
build_indexed_dataset_LDFLAGS = $(ADIOSLIB_LDFLAGS)

CLEANFILES = *.bp
CC=$(MPICC)
#EXTRA_DIST = test_config_flex.xml 

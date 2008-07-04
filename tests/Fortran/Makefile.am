INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS=adios_test_f

#adios_test_f_SOURCES = adios_test_f.F90 adios_module.F90
adios_test_f_SOURCES = adios_test_f.F90
adios_test_f_LDFLAGS = $(MXML_LDFLAGS)
adios_test_f_LDADD = $(top_builddir)/src/libadios.a
adios_test_f_LDADD += $(MXML_LIBS)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) "$()" \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

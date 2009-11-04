INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adiosf_global

adiosf_global_SOURCES = adiosf_global.F90
adiosf_global_LDADD = $(top_builddir)/src/libadios.a 
adiosf_global_LDFLAGS = $(MXML_LDFLAGS)
adiosf_global_LDADD += $(MXML_LIBS)
if HAVE_PHDF5
adiosf_global_LDFLAGS += $(PHDF5_LDFLAGS)
adiosf_global_LDADD +=  $(top_builddir)/src/libadios_phdf5.a $(PHDF5_LIBS)
endif
CLEANFILES = *.bp
CC=$(MPICC)

FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

EXTRA_DIST = adiosf_global.xml gread_temperature.fh gwrite_temperature.fh job.pbs

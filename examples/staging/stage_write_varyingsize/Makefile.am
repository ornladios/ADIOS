AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = stage_write 

stage_write_SOURCES =  stage_write.c utils.c decompose_block.c
stage_write_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
stage_write_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)


if BUILD_FORTRAN
    noinst_PROGRAMS += genarray_varying 
    genarray_varying_SOURCES = genarray_varying.F90
    genarray_varying_LDADD = $(top_builddir)/src/libadiosf.a $(ADIOSLIB_LDADD)
    genarray_varying_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
endif

CC=$(MPICC)
FC=$(MPIFC)

FCLINK = \
   $(LIBTOOL) --mode=link --tag F77 $(FC) \
   $(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

CLEANFILES = *.bp *.o[0-9]* *.mod *.fh core.* log.* dataspaces.conf conf srv.lck
EXTRA_DIST = utils.h decompose.h decompose_single.c test_decompose.c


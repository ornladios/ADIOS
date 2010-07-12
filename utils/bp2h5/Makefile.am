INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2h5

bp2h5_SOURCES = bp2h5.c hw-utils.c

bp2h5_CPPFLAGS = -DH5_USE_16_API $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSLIB_SEQ_CPPFLAGS) $(ADIOSLIB_SEQ_CFLAGS)
bp2h5_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS) -static-libtool-libs
bp2h5_LDADD = $(top_builddir)/src/libadios_internal_nompi.a
bp2h5_LDADD += $(ADIOSLIB_SEQ_LDADD)

EXTRA_DIST = hw-utils.h
if USE_PARALLEL_COMPILER
    CC=$(MPICC)
    CXX=$(MPICXX)
endif


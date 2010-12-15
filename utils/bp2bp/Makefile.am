INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2bp

bp2bp_SOURCES = bp2bp.c
bp2bp_CPPFLAGS = -DH5_USE_16_API -O0 -g $(ADIOSLIB_EXTRA_CPPFLAGS) $(ADIOSREADLIB_CPPFLAGS) $(ADIOSREADLIB_CFLAGS )
bp2bp_LDFLAGS = $(ADIOSREADLIB_LDFLAGS) 
bp2bp_LDADD = $(top_builddir)/src/libadiosread.a  $(top_builddir)/src/libadios.a

bp2bp_LDADD += $(ADIOSREADLIB_LDADD)

#CC=$(MPICC)
#CXX=$(MPICXX)



# Add your codes that you want to build at 'make check'
# Do not build Fortran codes if Fortran is not supported.
# Do not build a code using the write API if only the read API is available
# Steps:
# 1. add code name below to check_readonly_C/check_C/check_readonly_Fortran/check_Fortran targets
# 2. add code name to TESTS if it is a sequential test 
# 3. define build rules for your code
# 4. add files to CLEANFILES that should be deleted at 'make clean'
# 5. add to EXTRA_DIST any non-source files that should go with the distribution
check_readonly_C=
if BUILD_WRITE
check_C=local_array_time \
        adios_amr_write \
        adios_amr_write_2vars \
        adios_staged_read \
        adios_staged_read_v2 \
        adios_staged_read_2vars \
        write_read \
        write_alternate \
        two_groups \
        many_vars \
        selections \
        connect_to_space_subset
endif

if BUILD_FORTRAN
check_readonly_Fortran=
if BUILD_WRITE
check_Fortran=posix_method
endif
endif

check_PROGRAMS=$(check_C) $(check_Fortran) $(check_readonly_C) $(check_readonly_Fortran)

# Sequential tests to be executed by "make check" can be defined right here in TESTS
# Parallel tests should be defined in ../tests/*.sh and executed by ../test.sh
TESTS= 

INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

#############################
# Build rules for each test #
#############################

# posix_method
posix_method_SOURCES=posix_method.F90
posix_method_LDADD = $(top_builddir)/src/libadiosf.a $(ADIOSLIB_LDADD)
posix_method_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
posix_method.o: gwrite_posix_method.fh 
gwrite_posix_method.fh: posix_method.xml
	$(top_srcdir)/utils/gpp/gpp.py posix_method.xml

local_array_time_SOURCES=local_array_time.c
local_array_time_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
local_array_time_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)

adios_amr_write_SOURCES=adios_amr_write.c
adios_amr_write_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
adios_amr_write_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
adios_amr_write.o: adios_amr_write.c

adios_amr_write_2vars_SOURCES=adios_amr_write_2vars.c
adios_amr_write_2vars_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
adios_amr_write_2vars_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
adios_amr_write_2vars.o: adios_amr_write_2vars.c

adios_staged_read_SOURCES=adios_staged_read.c
adios_staged_read_LDADD = $(top_builddir)/src/libadiosread.a $(ADIOSREADLIB_LDADD)
adios_staged_read_LDFLAGS = $(AM_LDFLAGS) $(ADIOSREADLIB_LDFLAGS)
adios_staged_read.o: adios_staged_read.c

adios_staged_read_v2_SOURCES=adios_staged_read_v2.c
adios_staged_read_v2_LDADD = $(top_builddir)/src/libadiosread.a $(ADIOSREADLIB_LDADD)
adios_staged_read_v2_LDFLAGS = $(AM_LDFLAGS) $(ADIOSREADLIB_LDFLAGS)
adios_staged_read_v2.o: adios_staged_read_v2.c

adios_staged_read_2vars_SOURCES=adios_staged_read_2vars.c
adios_staged_read_2vars_LDADD = $(top_builddir)/src/libadiosread.a $(ADIOSREADLIB_LDADD)
adios_staged_read_2vars_LDFLAGS = $(AM_LDFLAGS) $(ADIOSREADLIB_LDFLAGS)
adios_staged_read_2vars.o: adios_staged_read_2vars.c

write_read_SOURCES=write_read.c
write_read_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
write_read_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
write_read.o: write_read.c

write_alternate_SOURCES=write_alternate.c
write_alternate_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
write_alternate_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
write_alternate.o: write_alternate.c

two_groups_SOURCES=two_groups.c
two_groups_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
two_groups_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
two_groups.o: two_groups.c

many_vars_SOURCES=many_vars.c
many_vars_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
many_vars_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
many_vars.o: many_vars.c

selections_SOURCES=selections.c
selections_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
selections_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
selections.o: selections.c

connect_to_space_subset_SOURCES=connect_to_space_subset.c
connect_to_space_subset_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
connect_to_space_subset_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_LDFLAGS)
connect_to_space_subset.o: connect_to_space_subset.c

EXTRA_DIST = adios_amr_write.xml adios_amr_write_2vars.xml \
             posix_method.xml local_array_time.xml  \
             write_alternate.xml write_read.xml

CLEANFILES = gwrite_posix_method.fh gread_posix_method.fh 

CC=$(MPICC)
FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

# This directory is for sequential tests of src/ functions only.
# Add your codes that you want to build at 'make check'
# Do not build Fortran codes if Fortran is not supported.
# Do not build a code using the write API if only the read API is available
# Steps:
# 1. add C code name to test_C
# 2. add Fortran code name to test_F
# 3. define build rules for your code
# 4. add files to CLEANFILES that should be deleted at 'make clean'
# 5. add to EXTRA_DIST any non-source files that should go with the distribution

test_C = hashtest copy_subvolume text_to_pairstruct test_strutil points_1DtoND trim_spaces

if BUILD_WRITE
    test_C += transforms_specparse group_free_test query_minmax read_points_2d read_points_3d
endif

if BUILD_FORTRAN
    test_F = selection_api
if BUILD_WRITE
#    test_F += 
endif
endif

check_PROGRAMS=$(test_C) $(test_F)

# Sequential tests to be executed by "make check" can be defined right here in TESTS
TESTS=$(test_C) $(test_F)

AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public  -I$(top_builddir)/src -I$(top_srcdir)/src/public

#all-local:
#	test "$(srcdir)" = "$(builddir)" || cp -p $(srcdir)/*.xml $(builddir)

#############################
# Build rules for each test #
#############################

#
# C Tests 
#
hashtest_SOURCES=hashtest.c
hashtest_LDADD = $(top_builddir)/src/core/libadios_a-qhashtbl.o
hashtest_LDFLAGS = $(AM_LDFLAGS) 
hashtest_CPPFLAGS = -I$(top_srcdir)/src
hashtest.o: hashtest.c

copy_subvolume_SOURCES=copy_subvolume.c
copy_subvolume_LDADD = $(top_builddir)/src/libadiosread_nompi.a $(ADIOSREADLIB_SEQ_LDADD)
copy_subvolume_LDFLAGS = $(AM_LDFLAGS) $(ADIOSREADLIB_SEQ_LDFLAGS)
copy_subvolume_CPPFLAGS = -I$(top_srcdir)/src $(ADIOSREADLIB_SEQ_CPPFLAGS) -I$(top_builddir)/src/public
copy_subvolume.o: copy_subvolume.c

text_to_pairstruct_SOURCES=text_to_pairstruct.c
text_to_pairstruct_LDADD = $(top_builddir)/src/core/libadios_a-strutil.o
text_to_pairstruct_LDFLAGS = $(AM_LDFLAGS) 
text_to_pairstruct_CPPFLAGS = -I$(top_srcdir)/src
text_to_pairstruct.o: text_to_pairstruct.c

test_strutil_SOURCES=test_strutil.c
test_strutil_LDADD = $(top_builddir)/src/core/libadios_a-strutil.o
test_strutil_LDFLAGS = $(AM_LDFLAGS) 
test_strutil_CPPFLAGS = -I$(top_srcdir)/src
test_strutil.o: test_strutil.c

points_1DtoND_SOURCES=points_1DtoND.c
points_1DtoND_LDADD = $(top_builddir)/src/libadiosread_nompi.a $(ADIOSREADLIB_SEQ_LDADD)
points_1DtoND_LDFLAGS = $(AM_LDFLAGS) $(ADIOSREADLIB_SEQ_LDFLAGS)
points_1DtoND_CPPFLAGS = -I$(top_srcdir)/src $(ADIOSREADLIB_SEQ_CPPFLAGS) -I$(top_builddir)/src/public
points_1DtoND.o: points_1DtoND.c

trim_spaces_SOURCES=trim_spaces.c
trim_spaces_LDADD = $(top_builddir)/src/libadiosread_nompi.a $(ADIOSREADLIB_SEQ_LDADD)
trim_spaces_LDFLAGS = $(AM_LDFLAGS) $(ADIOSREADLIB_SEQ_LDFLAGS)
trim_spaces_CPPFLAGS = -I$(top_srcdir)/src $(ADIOSREADLIB_SEQ_CPPFLAGS) -I$(top_builddir)/src/public
trim_spaces.o: trim_spaces.c

#
# C Tests built only with write-enabled
#
transforms_specparse_SOURCES=transforms_specparse.c
transforms_specparse_LDADD = $(top_builddir)/src/libadios_nompi.a $(ADIOSLIB_SEQ_LDADD)
transforms_specparse_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_SEQ_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
transforms_specparse_CPPFLAGS = -I$(top_srcdir)/src $(ADIOSLIB_SEQ_CPPFLAGS)
transforms_specparse.o: transforms_specparse.c

group_free_test_SOURCES=group_free_test.c
group_free_test_LDADD = $(top_builddir)/src/libadios_nompi.a $(ADIOSLIB_SEQ_LDADD)
group_free_test_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_SEQ_LDFLAGS)
group_free_test_CPPFLAGS = -I$(top_srcdir)/src $(ADIOSLIB_SEQ_CPPFLAGS) -I$(top_builddir)/src/public
group_free_test.o: group_free_test.c

query_minmax_SOURCES=query_minmax.c
query_minmax_LDADD = $(top_builddir)/src/libadios_nompi.a $(ADIOSLIB_SEQ_LDADD)
query_minmax_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_SEQ_LDFLAGS)
query_minmax_CPPFLAGS = -I$(top_srcdir)/src $(ADIOSLIB_SEQ_CPPFLAGS) -I$(top_builddir)/src/public
query_minmax.o: query_minmax.c

read_points_2d_SOURCES=read_points_2d.c
read_points_2d_LDADD = $(top_builddir)/src/libadios_nompi.a $(ADIOSLIB_SEQ_LDADD)
read_points_2d_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_SEQ_LDFLAGS)
read_points_2d_CPPFLAGS = -I$(top_srcdir)/src $(ADIOSLIB_SEQ_CPPFLAGS) -I$(top_builddir)/src/public
read_points_2d.o: read_points_2d.c

read_points_3d_SOURCES=read_points_3d.c
read_points_3d_LDADD = $(top_builddir)/src/libadios_nompi.a $(ADIOSLIB_SEQ_LDADD)
read_points_3d_LDFLAGS = $(AM_LDFLAGS) $(ADIOSLIB_SEQ_LDFLAGS)
read_points_3d_CPPFLAGS = -I$(top_srcdir)/src $(ADIOSLIB_SEQ_CPPFLAGS) -I$(top_builddir)/src/public
read_points_3d.o: read_points_3d.c

#
# FORTRAN Tests
#
selection_api_SOURCES=selection_api.F90
selection_api_LDADD = $(top_builddir)/src/libadiosreadf_nompi.a $(ADIOSREADLIB_SEQ_LDADD)
selection_api_LDFLAGS = $(AM_LDFLAGS) $(ADIOSREADLIB_SEQ_LDFLAGS)



EXTRA_DIST = 

CLEANFILES = 

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

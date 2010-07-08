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
check_C= 
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


EXTRA_DIST = posix_method.xml

CLEANFILES = gwrite_posix_method.fh gread_posix_method.fh

FC=$(MPIFC)

FCLINK = \
	$(LIBTOOL) --mode=link --tag F77 $(FC) \
	$(AM_FCFLAGS) $(FCFLAGS) $(LDFLAGS) -o $@

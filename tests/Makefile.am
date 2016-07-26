SUBDIRS=test_src
SUITEDIR=suite
if BUILD_WRITE
if  HAVE_MPI
    SUBDIRS += C bp_read 
if      BUILD_FORTRAN
        SUBDIRS += Fortran genarray
endif   BUILD_FORTRAN
    SUBDIRS += ${SUITEDIR}
endif HAVE_MPI
endif BUILD_WRITE

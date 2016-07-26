SUBDIRS=test_src
SUITEDIR=suite
if BUILD_WRITE
if  BUILD_MPI
    SUBDIRS += C bp_read 
if      BUILD_FORTRAN
        SUBDIRS += Fortran genarray
endif   BUILD_FORTRAN
    SUBDIRS += ${SUITEDIR}
endif BUILD_MPI
endif BUILD_WRITE

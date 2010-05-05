
AC_INIT(configure.ac)
AC_CONFIG_AUX_DIR(config)
AC_REVISION($Revision$)

AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE(adios,1.2)

AC_LANG(C)
AC_C_STRINGIZE
dnl generates HAVE_STRINGIZE

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX

AC_LANG(C)
ACX_MPI([], [AC_MSG_ERROR([couldn't find mpi library for C])])

AC_LANG(C++)
ACX_MPI([], [AC_MSG_ERROR([couldn't find mpi library for C++])])

AC_ARG_ENABLE(fortran, 
    [AS_HELP_STRING([--disable-fortran],[disable Fortran test programs in ADIOS.])])
AM_CONDITIONAL([BUILD_FORTRAN], [test "x$enable_fortran" != "xno"])

if test "x$enable_fortran" != "xno"; then
    AC_PROG_FC
    AC_PROG_FC_C_O
    AC_FC_WRAPPERS
    dnl AC_PROG_F77
    dnl AC_PROG_F77_C_O
    AC_LANG(Fortran)
    ACX_MPI([], [AC_MSG_ERROR([couldn't find mpi library for Fortran])])
fi

AC_LANG(C)
AC_PROG_LIBTOOL
AC_PROG_MXML
AC_PROG_GENGETOPT
AC_INFINIBAND
AC_PORTALS
AC_HDF5
AC_PHDF5
AC_NETCDF
AC_NSSI
AC_PROG_DART
AC_DMALLOC
AC_DATATAP

AC_LANG(C)


AC_ARG_ENABLE(research-transports, 
    [AS_HELP_STRING([--enable-research-transports], 
        [Enable building the various research transports. Default is not to build these as they are not part of the base distribution.])])
AM_CONDITIONAL([RESEARCH_TRANSPORTS], [test "x$enable_research_transports" == "xyes"])

AC_LANG(C)

AM_CONDITIONAL([HAVE_DATATAP], [test x$datatap != xdisable])

AC_ARG_ENABLE(write, 
    [AS_HELP_STRING([--disable-write],[disable building the write methods in ADIOS. Must be used on MacOSX])])
AM_CONDITIONAL([BUILD_WRITE], [test "x$enable_write" != "xno"])

dnl CC="$MPICC"
CXX="$MPICXX"
FC="$MPIFC"
CPPFLAGS="$CPPFLAGS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_GNU_SOURCE"

LIBS="$MPILIBS $LIBS"

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile
                 tests/C/Makefile
                 tests/Fortran/Makefile
                 tests/genarray/Makefile
                 tests/bp_read/Makefile
                 tests/suite/Makefile
                 tests/suite/programs/Makefile
                 utils/Makefile
                 utils/adios_lint/Makefile
                 utils/gpp/Makefile
                 utils/bpdump/Makefile
                 utils/bp2h5/Makefile
                 utils/bp2ncd/Makefile
                 utils/bp2ascii/Makefile
                 utils/bpsplit/Makefile
                 utils/bpls/Makefile
                 examples/Makefile
                 examples/C/Makefile
                 examples/C/scalars/Makefile
                 examples/C/arrays/Makefile
                 examples/C/attributes/Makefile
                 examples/C/manual/Makefile
                 examples/C/global-array/Makefile
                 examples/C/global-array-time/Makefile
                 examples/C/read_all/Makefile
                 examples/Fortran/Makefile
                 examples/Fortran/scalars/Makefile
                 examples/Fortran/arrays/Makefile
                 examples/Fortran/global-array/Makefile
                 examples/Fortran/global-array-time/Makefile
                ])

AC_OUTPUT

echo
echo "*********************************"
echo $PACKAGE configured as follows:
echo "*********************************"
echo

echo
echo "  - General Compile FLAGS";
echo "      - FC = $FC";
echo "      - FCFLAGS = $FCFLAGS";
echo
echo "      - CC = $CC";
echo "      - CPPFLAGS = $CPPFLAGS";
echo "      - CFLAGS = $CFLAGS";
echo "      - LDLAGS = $LDFLAGS";
echo "      - LIBS = $LIBS";
echo
echo "      - MPICC = $MPICC";
echo "      - MPILIBS = $MPILIBS";
echo "      - MPIFC = $MPIFC";
echo

if test "x$enable_research_transports" == "xyes"; then 
        echo "  - Research Transports Included";

else
        echo "  - No research transport methods";
fi
echo

if test "x${datatap}" != "xdisable"; then
        echo "  - DataTap with interconnect: $datatap";
else	
        echo "  - No DataTap transport";
fi

if test -z "${HAVE_DART_TRUE}"; then
        echo "  - DART transport method";
        echo "      - DART_CFLAGS = $DART_CFLAGS";
        echo "      - DART_CPPFLAGS = $DART_CPPFLAGS";
        echo "      - DART_LDFLAGS = $DART_LDFLAGS";
        echo "      - DART_LIBS = $DART_LIBS";
else
	echo "  - No DART to build DART transport method"
fi
echo

dnl if test -z "${HAVE_PORTALS_TRUE}"; then
dnl        echo "  - PORTALS";
dnl        echo "      - PORTALS_CFLAGS = $PORTALS_CFLAGS";
dnl        echo "      - PORTALS_CPPFLAGS = $PORTALS_CPPFLAGS";
dnl        echo "      - PORTALS_LDFLAGS = $PORTALS_LDFLAGS";
dnl        echo "      - PORTALS_LIBS = $PORTALS_LIBS";
dnl else
dnl        echo " - No PORTALS"
dnl fi
dnl
dnl echo

if test -z "${HAVE_MXML_TRUE}"; then
        echo "  - MXML";
        echo "      - MXML_CFLAGS = $MXML_CFLAGS";
        echo "      - MXML_CPPFLAGS = $MXML_CPPFLAGS";
        echo "      - MXML_LDFLAGS = $MXML_LDFLAGS";
        echo "      - MXML_LIBS = $MXML_LIBS";
else
        echo "  - No MXML"
fi

echo

if test -z "${HAVE_GENGETOPT_TRUE}"; then
        echo "  - GENGETOPT";
        echo "      - GENGETOPT_CFLAGS = $GENGETOPT_CFLAGS";
        echo "      - GENGETOPT_CPPFLAGS = $GENGETOPT_CPPFLAGS";
        echo "      - GENGETOPT_LDFLAGS = $GENGETOPT_LDFLAGS";
        echo "      - GENGETOPT_LIBS = $GENGETOPT_LIBS";
else
        echo "  - No Gengetopt"
fi

echo

if test -z "${HAVE_MPI_TRUE}"; then
        echo "  - MPI";
        echo "      - MPI_CFLAGS = $MPI_CFLAGS";
        echo "      - MPI_CPPFLAGS = $MPI_CPPFLAGS";
        echo "      - MPI_LDFLAGS = $MPI_LDFLAGS";
        echo "      - MPI_LIBS = $MPI_LIBS";
else
        echo "  - No MPI"
fi

echo

if test -z "${HAVE_PHDF5_TRUE}"; then
        echo "  - PHDF5";
        echo "      - PHDF5_CFLAGS = $PHDF5_CFLAGS";
        echo "      - PHDF5_CPPFLAGS = $PHDF5_CPPFLAGS";
        echo "      - PHDF5_LDFLAGS = $PHDF5_LDFLAGS";
        echo "      - PHDF5_LIBS = $PHDF5_LIBS";
else
	echo "  - No PHDF5"
fi

echo

if test -z "${HAVE_HDF5_TRUE}"; then
        echo "  - HDF5";
        echo "      - HDF5_CFLAGS = $HDF5_CFLAGS";
        echo "      - HDF5_CPPFLAGS = $HDF5_CPPFLAGS";
        echo "      - HDF5_LDFLAGS = $HDF5_LDFLAGS";
        echo "      - HDF5_LIBS = $HDF5_LIBS";
else
	echo "  - No HDF5"
fi

echo

if test -z "${HAVE_NETCDF_TRUE}"; then
        echo "  - NetCDF";
        echo "      - NETCDF_CFLAGS = $NETCDF_CFLAGS";
        echo "      - NETCDF_CPPFLAGS = $NETCDF_CPPFLAGS";
        echo "      - NETCDF_LDFLAGS = $NETCDF_LDFLAGS";
        echo "      - NETCDF_LIBS = $NETCDF_LIBS";
else
	echo "  - No NETCDF"
fi

echo

if test -z "${HAVE_PORTALS_TRUE}"; then
        echo "  - PORTALS";
        echo "      - PORTALS_CFLAGS = $PORTALS_CFLAGS";
        echo "      - PORTALS_CPPFLAGS = $PORTALS_CPPFLAGS";
        echo "      - PORTALS_LDFLAGS = $PORTALS_LDFLAGS";
        echo "      - PORTALS_LIBS = $PORTALS_LIBS";
else
	echo "  - No Portals"
fi

echo

if test -z "${HAVE_INFINIBAND_TRUE}"; then
        echo "  - INFINIBAND";
        echo "      - INFINIBAND_CFLAGS = $INFINIBAND_CFLAGS";
        echo "      - INFINIBAND_CPPFLAGS = $INFINIBAND_CPPFLAGS";
        echo "      - INFINIBAND_LDFLAGS = $INFINIBAND_LDFLAGS";
        echo "      - INFINIBAND_LIBS = $INFINIBAND_LIBS";
else
	echo "  - No IB"
fi

echo

if test -z "${HAVE_NSSI_TRUE}"; then
        echo "  - nessie";
        echo "      - NSSI_CFLAGS = $NSSI_CFLAGS";
        echo "      - NSSI_CPPFLAGS = $NSSI_CPPFLAGS";
        echo "      - NSSI_LDFLAGS = $NSSI_LDFLAGS";
        echo "      - NSSI_LIBS = $NSSI_LIBS";
else
	echo "  - No NSSI"

fi

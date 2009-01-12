
AC_INIT(configure.ac)
AC_CONFIG_AUX_DIR(config)
AC_REVISION($Revision$)

AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE(ADIOS,0.8)

AC_LANG(C)
AC_C_STRINGIZE
dnl generates HAVE_STRINGIZE

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX

AC_PROG_FC
AC_PROG_FC_C_O

dnl AC_PROG_F77
dnl AC_PROG_F77_C_O

AC_PROG_LIBTOOL

AC_PROG_MXML
AC_PROG_GENGETOPT
AC_HDF5
AC_PHDF5
AC_NETCDF

AC_LANG(C)
ACX_MPI([], [AC_MSG_ERROR([couldn't find mpi library for C])])

AC_LANG(Fortran)
ACX_MPI([], [AC_MSG_ERROR([couldn't find mpi library for Fortran])])

AC_LANG(C)

CC="$MPICC"
CXX="$MPICXX"
FC="$MPIFC"

LIBS="$MPILIBS $LIBS"

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile
                 tests/C/Makefile
                 tests/Fortran/Makefile
                 utils/Makefile
                 utils/adios_lint/Makefile
                 utils/gpp/Makefile
                 utils/bpdump/Makefile
                 utils/bp2h5/Makefile
                 utils/bp2ncd/Makefile
                 utils/bp2ascii/Makefile
                 examples-code/Makefile
                 examples-code/global-array/Makefile
                 examples-code/global-array-time/Makefile
                 examples-code/time-index/Makefile
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

if test -z "${HAVE_PORTALS_TRUE}"; then
        echo "  - PORTALS";
        echo "      - PORTALS_CFLAGS = $PORTALS_CFLAGS";
        echo "      - PORTALS_CPPFLAGS = $PORTALS_CPPFLAGS";
        echo "      - PORTALS_LDFLAGS = $PORTALS_LDFLAGS";
        echo "      - PORTALS_LIBS = $PORTALS_LIBS";
else
echo "No Portals"
fi

echo

if test -z "${HAVE_MXML_TRUE}"; then
        echo "  - MXML";
        echo "      - MXML_CFLAGS = $MXML_CFLAGS";
        echo "      - MXML_CPPFLAGS = $MXML_CPPFLAGS";
        echo "      - MXML_LDFLAGS = $MXML_LDFLAGS";
        echo "      - MXML_LIBS = $MXML_LIBS";
else
echo "No MXML"
fi

echo

if test -z "${HAVE_GENGETOPT_TRUE}"; then
        echo "  - GENGETOPT";
        echo "      - GENGETOPT_CFLAGS = $GENGETOPT_CFLAGS";
        echo "      - GENGETOPT_CPPFLAGS = $GENGETOPT_CPPFLAGS";
        echo "      - GENGETOPT_LDFLAGS = $GENGETOPT_LDFLAGS";
        echo "      - GENGETOPT_LIBS = $GENGETOPT_LIBS";
else
echo "No Gengetopt"
fi

echo

if test -z "${HAVE_MPI_TRUE}"; then
        echo "  - MPI";
        echo "      - MPI_CFLAGS = $MPI_CFLAGS";
        echo "      - MPI_CPPFLAGS = $MPI_CPPFLAGS";
        echo "      - MPI_LDFLAGS = $MPI_LDFLAGS";
        echo "      - MPI_LIBS = $MPI_LIBS";
else
echo "No MPI"
fi

echo

if test -z "${HAVE_PHDF5_TRUE}"; then
        echo "  - PHDF5";
        echo "      - PHDF5_CFLAGS = $PHDF5_CFLAGS";
        echo "      - PHDF5_CPPFLAGS = $PHDF5_CPPFLAGS";
        echo "      - PHDF5_LDFLAGS = $PHDF5_LDFLAGS";
        echo "      - PHDF5_LIBS = $PHDF5_LIBS";
else
	echo "No PHDF5"
fi

echo

if test -z "${HAVE_HDF5_TRUE}"; then
        echo "  - HDF5";
        echo "      - HDF5_CFLAGS = $HDF5_CFLAGS";
        echo "      - HDF5_CPPFLAGS = $HDF5_CPPFLAGS";
        echo "      - HDF5_LDFLAGS = $HDF5_LDFLAGS";
        echo "      - HDF5_LIBS = $HDF5_LIBS";
else
	echo "No HDF5"
fi

echo

if test -z "${HAVE_NETCDF_TRUE}"; then
        echo "  - NetCDF";
        echo "      - NETCDF_CFLAGS = $NETCDF_CFLAGS";
        echo "      - NETCDF_CPPFLAGS = $NETCDF_CPPFLAGS";
        echo "      - NETCDF_LDFLAGS = $NETCDF_LDFLAGS";
        echo "      - NETCDF_LIBS = $NETCDF_LIBS";
else
	echo "No NETCDF"
fi



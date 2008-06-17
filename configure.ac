
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
                ])

dnl add these in as they are converted: [utils/bp2ascii/Makefile utils/bp2h5/Makefile utils/bp2ncd/Makefile utils/gwrite/bpdump utils/gwrite/Makefile])

AC_OUTPUT

echo
echo "*********************************"
echo $PACKAGE configured as follows:
echo "*********************************"
echo

if test -z "${HAVE_PORTALS_TRUE}"; then
echo "Has Portals"
else
echo "No Portals"
fi

echo

if test -z "${HAVE_MXML_TRUE}"; then
echo "Has MXML"
else
echo "No MXML"
fi

echo

if test -z "${HAVE_GENGETOPT_TRUE}"; then
echo "has Gengetopt"
else
echo "No Gengetopt"
fi

echo

if test -z "${HAVE_MPI_TRUE}"; then
echo "has MPI"
else
echo "No MPI"
fi


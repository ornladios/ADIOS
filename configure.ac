
AC_INIT(configure.ac)
AC_CONFIG_AUX_DIR(config)
AC_REVISION($Revision$)

AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE(adios,1.2.2)

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
    F77=$FC
    dnl AC_PROG_F77
    dnl AC_PROG_F77_C_O
    AC_LANG(Fortran)
    ACX_MPI([], [AC_MSG_ERROR([couldn't find mpi library for Fortran])])
fi

AC_LANG(C)

dnl Generate libtool script and then do a nasty thing.
dnl PGI compiler (behind MPICC) does not like --rpath.
dnl The only way I found to avoid the problem is to remove the chance of using it
dnl inside the libtool script. Similar to what HDF5 configure does.
dnl see: http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/configure.in
dnl look the snippet after LT_INIT there.
AC_PROG_LIBTOOL
#LT_INIT([disable-shared])
dnl echo 'fixing hardcode_libdir_flag_spec in' $ofile
dnl ed - $ofile <<EOF 2> /dev/null
dnl g/^wl=".*"/s//wl="-Wl,"/
dnl w
dnl q
dnl EOF
dnl g/^hardcode_libdir_flag_spec=".*"/s//hardcode_libdir_flag_spec=""/
dnl g/^hardcode_libdir_flag_spec_ld=".*"/s//hardcode_libdir_flag_spec_ld=""/

AC_PROG_MXML
dnl AC_PROG_GENGETOPT
AC_INFINIBAND
AC_PORTALS
AC_HDF5
AC_PHDF5
AC_NETCDF
AC_NC4PAR
AC_NSSI
AC_PROG_DART
AC_PROG_DIMES
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
dnl CXX="$MPICXX"
dnl FC="$MPIFC"

LIBS="$MPILIBS"

AM_CONDITIONAL(USE_PARALLEL_COMPILER,false)

dnl Use PHDF5 for bp2h5 if no sequential HDF5 was found
if test -z "${HAVE_PHDF5_TRUE}" -a -z "${HAVE_HDF5_FALSE}"; then
    HDF5_CPPFLAGS="${PHDF5_CPPFLAGS}"
    HDF5_CFLAGS="${PHDF5_CFLAGS}"
    HDF5_LDFLAGS="${PHDF5_LDFLAGS}"
    HDF5_LIBS="${PHDF5_LIBS}"
    AM_CONDITIONAL(HAVE_HDF5,true)
    AC_DEFINE(HAVE_HDF5,1,[Define if you have HDF5.])
    AM_CONDITIONAL(USE_PARALLEL_COMPILER,true)
    AC_MSG_NOTICE([Use PHDF5 for sequential HDF5 tools])
fi

dnl Use NC4PAR for bp2ncd if no sequential NETCDF was found
if test -z "${HAVE_NC4PAR_TRUE}" -a -z "${HAVE_NETCDF_FALSE}"; then
    NETCDF_CPPFLAGS="${NC4PAR_CPPFLAGS}"
    NETCDF_CFLAGS="${NC4PAR_CFLAGS}"
    NETCDF_LDFLAGS="${NC4PAR_LDFLAGS}"
    NETCDF_LIBS="${NC4PAR_LIBS}"
    AM_CONDITIONAL(HAVE_NETCDF,true)
    AC_DEFINE(HAVE_NETCDF,1,[Define if you have NETCDF.])
    AM_CONDITIONAL(USE_PARALLEL_COMPILER,true)
    AC_MSG_NOTICE([Use NC4PAR for sequential NETCDF tools])
fi


dnl Generate flags for adios lib creation which depends on the packages we found
dnl ADIOSLIB_* for compiling libadios.a and libadiosf.a, and linking codes with those
dnl ADIOSLIB_SEQ_* for compiling libadios_nompi.a and libadiosf_nompi.a with dummy MPI
dnl ADIOSLIB_INT_* for codes using libadios_internal_nompi.a (like bp2ncd, bp2h5, bpdump etc)
dnl ADIOSREADLIB_* for codes using libadiosread.a, libadiosreadf.a
dnl ADIOSREADLIB_SEQ_* for codes using libadiosread_nompi.a
dnl adios_config will contain flags for users' code linking based on these flags
dnl EXTRA flags are used for compiling adios stuff, which need not go into adios_config for user codes
ADIOSLIB_EXTRA_CPPFLAGS="-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_GNU_SOURCE"
ADIOSLIB_CPPFLAGS="${MXML_CPPFLAGS}"
ADIOSLIB_CFLAGS=
ADIOSLIB_LDFLAGS="${MXML_LDFLAGS}"
ADIOSLIB_LDADD="-lm ${MXML_LIBS}"
ADIOSLIB_SEQ_CPPFLAGS="-D_NOMPI ${MXML_CPPFLAGS}"
ADIOSLIB_SEQ_CFLAGS=
ADIOSLIB_SEQ_LDFLAGS="${MXML_LDFLAGS}"
ADIOSLIB_SEQ_LDADD="-lm ${MXML_LIBS}"
ADIOSLIB_INT_CPPFLAGS="-D_NOMPI ${MXML_CPPFLAGS}"
ADIOSLIB_INT_CFLAGS=
ADIOSLIB_INT_LDFLAGS="${MXML_LDFLAGS}"
ADIOSLIB_INT_LDADD="-lm ${MXML_LIBS}"
ADIOSREADLIB_CPPFLAGS=
ADIOSREADLIB_CFLAGS=
ADIOSREADLIB_LDFLAGS=
ADIOSREADLIB_LDADD="-lm"
ADIOSREADLIB_SEQ_CPPFLAGS="-D_NOMPI"
ADIOSREADLIB_SEQ_CFLAGS=
ADIOSREADLIB_SEQ_LDFLAGS=
ADIOSREADLIB_SEQ_LDADD="-lm"
if test "x${datatap}" != "xdisable"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} ${DT_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS} ${DT_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${DT_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${DT_LIBS}"
    ADIOSREADLIB_CPPFLAGS="${ADIOSREADLIB_CPPFLAGS} ${DT_CPPFLAGS}"
    ADIOSREADLIB_CFLAGS="${ADIOSREADLIB_CFLAGS} ${DT_CFLAGS}"
    ADIOSREADLIB_LDFLAGS="${ADIOSREADLIB_LDFLAGS} ${DT_LDFLAGS}"
    ADIOSREADLIB_LDADD="${ADIOSREADLIB_LDADD} ${DT_LIBS}"
    ADIOSREADLIB_SEQ_CPPFLAGS="${ADIOSREADLIB_SEQ_CPPFLAGS} ${DT_CPPFLAGS}"
    ADIOSREADLIB_SEQ_CFLAGS="${ADIOSREADLIB_SEQ_CFLAGS} ${DT_CFLAGS}"
    ADIOSREADLIB_SEQ_LDFLAGS="${ADIOSREADLIB_SEQ_LDFLAGS} ${DT_LDFLAGS}"
    ADIOSREADLIB_SEQ_LDADD="${ADIOSREADLIB_SEQ_LDADD} ${DT_LIBS}"
fi
if test -z "${HAVE_DART_TRUE}"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} ${DART_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${DART_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${DART_LIBS}"
    ADIOSLIB_SEQ_CPPFLAGS="${ADIOSLIB_SEQ_CPPFLAGS} ${DART_CPPFLAGS}"
    ADIOSLIB_SEQ_CFLAGS="${ADIOSLIB_SEQ_CFLAGS}"
    ADIOSLIB_SEQ_LDFLAGS="${ADIOSLIB_SEQ_LDFLAGS} ${DART_LDFLAGS}"
    ADIOSLIB_SEQ_LDADD="${ADIOSLIB_SEQ_LDADD} ${DART_LIBS}"
    ADIOSREADLIB_CPPFLAGS="${ADIOSREADLIB_CPPFLAGS} ${DART_CPPFLAGS}"
    ADIOSREADLIB_CFLAGS="${ADIOSREADLIB_CFLAGS} ${DART_CFLAGS}"
    ADIOSREADLIB_LDFLAGS="${ADIOSREADLIB_LDFLAGS} ${DART_LDFLAGS}"
    ADIOSREADLIB_LDADD="${ADIOSREADLIB_LDADD} ${DART_LIBS}"
    ADIOSREADLIB_SEQ_CPPFLAGS="${ADIOSREADLIB_SEQ_CPPFLAGS} ${DART_CPPFLAGS}"
    ADIOSREADLIB_SEQ_CFLAGS="${ADIOSREADLIB_SEQ_CFLAGS} ${DART_CFLAGS}"
    ADIOSREADLIB_SEQ_LDFLAGS="${ADIOSREADLIB_SEQ_LDFLAGS} ${DART_LDFLAGS}"
    ADIOSREADLIB_SEQ_LDADD="${ADIOSREADLIB_SEQ_LDADD} ${DART_LIBS}"
fi
if test -z "${HAVE_DIMES_TRUE}"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} ${DIMES_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${DIMES_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${DIMES_LIBS}"
    ADIOSLIB_SEQ_CPPFLAGS="${ADIOSLIB_SEQ_CPPFLAGS} ${DIMES_CPPFLAGS}"
    ADIOSLIB_SEQ_CFLAGS="${ADIOSLIB_SEQ_CFLAGS}"
    ADIOSLIB_SEQ_LDFLAGS="${ADIOSLIB_SEQ_LDFLAGS} ${DIMES_LDFLAGS}"
    ADIOSLIB_SEQ_LDADD="${ADIOSLIB_SEQ_LDADD} ${DIMES_LIBS}"
    ADIOSREADLIB_CPPFLAGS="${ADIOSREADLIB_CPPFLAGS} ${DIMES_CPPFLAGS}"
    ADIOSREADLIB_CFLAGS="${ADIOSREADLIB_CFLAGS} ${DIMES_CFLAGS}"
    ADIOSREADLIB_LDFLAGS="${ADIOSREADLIB_LDFLAGS} ${DIMES_LDFLAGS}"
    ADIOSREADLIB_LDADD="${ADIOSREADLIB_LDADD} ${DIMES_LIBS}"
    ADIOSREADLIB_SEQ_CPPFLAGS="${ADIOSREADLIB_SEQ_CPPFLAGS} ${DIMES_CPPFLAGS}"
    ADIOSREADLIB_SEQ_CFLAGS="${ADIOSREADLIB_SEQ_CFLAGS} ${DIMES_CFLAGS}"
    ADIOSREADLIB_SEQ_LDFLAGS="${ADIOSREADLIB_SEQ_LDFLAGS} ${DIMES_LDFLAGS}"
    ADIOSREADLIB_SEQ_LDADD="${ADIOSREADLIB_SEQ_LDADD} ${DIMES_LIBS}"
fi
if test -z "${HAVE_NSSI_TRUE}"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} ${NSSI_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS} ${NSSI_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${NSSI_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${NSSI_LIBS}"
    ADIOSREADLIB_CPPFLAGS="${ADIOSREADLIB_CPPFLAGS} ${MXML_CPPFLAGS} ${NSSI_CPPFLAGS}"
    ADIOSREADLIB_CFLAGS="${ADIOSREADLIB_CFLAGS} ${MXML_CPPFLAGS} ${NSSI_CFLAGS}"
    ADIOSREADLIB_LDFLAGS="${ADIOSREADLIB_LDFLAGS} ${MXML_LDFLAGS} ${NSSI_LDFLAGS}"
    ADIOSREADLIB_LDADD="${ADIOSREADLIB_LDADD} ${MXML_LIBS} ${NSSI_LIBS}"
    ADIOSREADLIB_SEQ_CPPFLAGS="${ADIOSREADLIB_SEQ_CPPFLAGS} ${MXML_CPPFLAGS} ${NSSI_CPPFLAGS}"
    ADIOSREADLIB_SEQ_CFLAGS="${ADIOSREADLIB_SEQ_CFLAGS} ${MXML_CPPFLAGS} ${NSSI_CFLAGS}"
    ADIOSREADLIB_SEQ_LDFLAGS="${ADIOSREADLIB_SEQ_LDFLAGS} ${MXML_LDFLAGS} ${NSSI_LDFLAGS}"
    ADIOSREADLIB_SEQ_LDADD="${ADIOSREADLIB_SEQ_LDADD} ${MXML_LIBS} ${NSSI_LIBS}"
fi
if test -z "${HAVE_NC4PAR_TRUE}"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} ${NC4PAR_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS} ${NC4PAR_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${NC4PAR_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${NC4PAR_LIBS}"
fi
if test -z "${HAVE_PHDF5_TRUE}"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} ${PHDF5_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS} ${PHDF5_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${PHDF5_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${PHDF5_LIBS}"
fi
if test -z "${HAVE_PORTALS_TRUE}"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} ${PORTALS_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS} ${PORTALS_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${PORTALS_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${PORTALS_LIBS}"
    ADIOSLIB_SEQ_CPPFLAGS="${ADIOSLIB_SEQ_CPPFLAGS} ${PORTALS_CPPFLAGS}"
    ADIOSLIB_SEQ_CFLAGS="${ADIOSLIB_SEQ_CFLAGS} ${PORTALS_CFLAGS}"
    ADIOSLIB_SEQ_LDFLAGS="${ADIOSLIB_SEQ_LDFLAGS} ${PORTALS_LDFLAGS}"
    ADIOSLIB_SEQ_LDADD="${ADIOSLIB_SEQ_LDADD} ${PORTALS_LIBS}"
fi
if test -z "${HAVE_INFINIBAND_TRUE}"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} ${INFINIBAND_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS} ${INFINIBAND_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${INFINIBAND_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${INFINIBAND_LIBS}"
    ADIOSLIB_SEQ_CPPFLAGS="${ADIOSLIB_SEQ_CPPFLAGS} ${INFINIBAND_CPPFLAGS}"
    ADIOSLIB_SEQ_CFLAGS="${ADIOSLIB_SEQ_CFLAGS} ${INFINIBAND_CFLAGS}"
    ADIOSLIB_SEQ_LDFLAGS="${ADIOSLIB_SEQ_LDFLAGS} ${INFINIBAND_LDFLAGS}"
    ADIOSLIB_SEQ_LDADD="${ADIOSLIB_SEQ_LDADD} ${INFINIBAND_LIBS}"
fi
if test -z "${HAVE_NETCDF_TRUE}"; then
    ADIOSLIB_INT_CPPFLAGS="${ADIOSLIB_INT_CPPFLAGS} ${NETCDF_CPPFLAGS}"
    ADIOSLIB_INT_CFLAGS="${ADIOSLIB_INT_CFLAGS} ${NETCDF_CFLAGS}"
    ADIOSLIB_INT_LDFLAGS="${ADIOSLIB_INT_LDFLAGS} ${NETCDF_LDFLAGS}"
    ADIOSLIB_INT_LDADD="${ADIOSLIB_INT_LDADD} ${NETCDF_LIBS}"
fi
if test -z "${HAVE_HDF5_TRUE}"; then
    ADIOSLIB_INT_CPPFLAGS="${ADIOSLIB_INT_CPPFLAGS} ${HDF5_CPPFLAGS}"
    ADIOSLIB_INT_CFLAGS="${ADIOSLIB_INT_CFLAGS} ${HDF5_CFLAGS}"
    ADIOSLIB_INT_LDFLAGS="${ADIOSLIB_INT_LDFLAGS} ${HDF5_LDFLAGS}"
    ADIOSLIB_INT_LDADD="${ADIOSLIB_INT_LDADD} ${HDF5_LIBS}"
fi
if test -z "${HAVE_DMALLOC_TRUE}"; then
    ADIOSLIB_CPPFLAGS="${ADIOSLIB_CPPFLAGS} -DDMALLOC ${DMALLOC_CPPFLAGS}"
    ADIOSLIB_CFLAGS="${ADIOSLIB_CFLAGS} ${DMALLOC_CFLAGS}"
    ADIOSLIB_LDFLAGS="${ADIOSLIB_LDFLAGS} ${DMALLOC_LDFLAGS}"
    ADIOSLIB_LDADD="${ADIOSLIB_LDADD} ${DMALLOC_LIBS}"
    ADIOSLIB_SEQ_CPPFLAGS="${ADIOSLIB_SEQ_CPPFLAGS} ${DMALLOC_CPPFLAGS}"
    ADIOSLIB_SEQ_CFLAGS="${ADIOSLIB_SEQ_CFLAGS} ${DMALLOC_CFLAGS}"
    ADIOSLIB_SEQ_LDFLAGS="${ADIOSLIB_SEQ_LDFLAGS} ${DMALLOC_LDFLAGS}"
    ADIOSLIB_SEQ_LDADD="${ADIOSLIB_SEQ_LDADD} ${DMALLOC_LIBS}"
    ADIOSLIB_INT_CPPFLAGS="${ADIOSLIB_INT_CPPFLAGS} ${DMALLOC_CPPFLAGS}"
    ADIOSLIB_INT_CFLAGS="${ADIOSLIB_INT_CFLAGS} ${DMALLOC_CFLAGS}"
    ADIOSLIB_INT_LDFLAGS="${ADIOSLIB_INT_LDFLAGS} ${DMALLOC_LDFLAGS}"
    ADIOSLIB_INT_LDADD="${ADIOSLIB_INT_LDADD} ${DMALLOC_LIBS}"
    ADIOSREADLIB_CPPFLAGS="${ADIOSREADLIB_CPPFLAGS} ${DMALLOC_CPPFLAGS}"
    ADIOSREADLIB_CFLAGS="${ADIOSREADLIB_CFLAGS} ${DMALLOC_CFLAGS}"
    ADIOSREADLIB_LDFLAGS="${ADIOSREADLIB_LDFLAGS} ${DMALLOC_LDFLAGS}"
    ADIOSREADLIB_LDADD="${ADIOSREADLIB_LDADD} ${DMALLOC_LIBS}"
    ADIOSREADLIB_SEQ_CPPFLAGS="${ADIOSREADLIB_SEQ_CPPFLAGS} ${DMALLOC_CPPFLAGS}"
    ADIOSREADLIB_SEQ_CFLAGS="${ADIOSREADLIB_SEQ_CFLAGS} ${DMALLOC_CFLAGS}"
    ADIOSREADLIB_SEQ_LDFLAGS="${ADIOSREADLIB_SEQ_LDFLAGS} ${DMALLOC_LDFLAGS}"
    ADIOSREADLIB_SEQ_LDADD="${ADIOSREADLIB_SEQ_LDADD} ${DMALLOC_LIBS}"
fi
AC_SUBST(ADIOSLIB_CPPFLAGS)
AC_SUBST(ADIOSLIB_CFLAGS)
AC_SUBST(ADIOSLIB_LDFLAGS)
AC_SUBST(ADIOSLIB_LDADD)
AC_SUBST(ADIOSLIB_SEQ_CPPFLAGS)
AC_SUBST(ADIOSLIB_SEQ_CFLAGS)
AC_SUBST(ADIOSLIB_SEQ_LDFLAGS)
AC_SUBST(ADIOSLIB_SEQ_LDADD)
AC_SUBST(ADIOSLIB_INT_CPPFLAGS)
AC_SUBST(ADIOSLIB_INT_CFLAGS)
AC_SUBST(ADIOSLIB_INT_LDFLAGS)
AC_SUBST(ADIOSLIB_INT_LDADD)
AC_SUBST(ADIOSREADLIB_CPPFLAGS)
AC_SUBST(ADIOSREADLIB_CFLAGS)
AC_SUBST(ADIOSREADLIB_LDFLAGS)
AC_SUBST(ADIOSREADLIB_LDADD)
AC_SUBST(ADIOSREADLIB_SEQ_CPPFLAGS)
AC_SUBST(ADIOSREADLIB_SEQ_CFLAGS)
AC_SUBST(ADIOSREADLIB_SEQ_LDFLAGS)
AC_SUBST(ADIOSREADLIB_SEQ_LDADD)
AC_SUBST(ADIOSLIB_EXTRA_CPPFLAGS)

if test -z "${USE_PARALLEL_COMPILER_TRUE}"; then
        AC_DEFINE(USE_PARALLEL_COMPILER, 1, [Need to use MPI compiler for sequential utils])
        echo "#### Defined USE_PARALLEL_COMPILER to = [${USE_PARALLEL_COMPILER}]"
fi

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
                 utils/bp2bp/Makefile
                 examples/Makefile
                 examples/C/Makefile
                 examples/C/scalars/Makefile
                 examples/C/arrays/Makefile
                 examples/C/attributes/Makefile
                 examples/C/manual/Makefile
                 examples/C/global-array/Makefile
                 examples/C/global-array-time/Makefile
                 examples/C/stat/Makefile
                 examples/C/read_all/Makefile
                 examples/Fortran/Makefile
                 examples/Fortran/scalars/Makefile
                 examples/Fortran/arrays/Makefile
                 examples/Fortran/global-array/Makefile
                 examples/Fortran/global-array-time/Makefile
                 adios_config.flags
                ])

AC_OUTPUT

# Get the version of the source we are configuring here
# Just to print it out at the end of configure
# see also Makefile.am:install-data-hook
if type svnversion &>/dev/null; then
    SVNVERSION=`svnversion`
else
    SVNVERSION="unknown"
fi

echo
echo "*********************************"
echo $PACKAGE configured as follows:
echo "*********************************"
echo
echo "  - General Compile FLAGS";
echo "      - CC  = $CC";
echo "      - CXX = $CXX";
echo "      - FC  = $FC";
echo "      - CPPFLAGS = $CPPFLAGS";
echo "      - CFLAGS   = $CFLAGS";
echo "      - CXXFLAGS = $CXXFLAGS";
echo "      - FCFLAGS  = $FCFLAGS";
echo "      - LDLAGS   = $LDFLAGS";
echo "      - LIBS     = $LIBS";
echo
echo "      - MPICC   = $MPICC";
echo "      - MPICXX  = $MPICXX";
echo "      - MPIFC   = $MPIFC";
echo "      - MPILIBS = $MPILIBS";
echo

if test "x$enable_research_transports" == "xyes"; then
    echo "  - Research Transports Included";

else
    echo "  - No research transport methods";
fi

if test -z "${HAVE_MXML_TRUE}"; then
    echo "  - MXML";
    echo "      - MXML_CFLAGS = $MXML_CFLAGS";
    echo "      - MXML_CPPFLAGS = $MXML_CPPFLAGS";
    echo "      - MXML_LDFLAGS = $MXML_LDFLAGS";
    echo "      - MXML_LIBS = $MXML_LIBS";
    echo
else
    echo "  - No MXML"
fi

if test -z "${HAVE_PORTALS_TRUE}"; then
    echo "  - PORTALS";
    echo "      - PORTALS_CFLAGS = $PORTALS_CFLAGS";
    echo "      - PORTALS_CPPFLAGS = $PORTALS_CPPFLAGS";
    echo "      - PORTALS_LDFLAGS = $PORTALS_LDFLAGS";
    echo "      - PORTALS_LIBS = $PORTALS_LIBS";
    echo
else
    echo "  - No Portals"
fi

if test -z "${HAVE_INFINIBAND_TRUE}"; then
    echo "  - INFINIBAND";
    echo "      - INFINIBAND_CFLAGS = $INFINIBAND_CFLAGS";
    echo "      - INFINIBAND_CPPFLAGS = $INFINIBAND_CPPFLAGS";
    echo "      - INFINIBAND_LDFLAGS = $INFINIBAND_LDFLAGS";
    echo "      - INFINIBAND_LIBS = $INFINIBAND_LIBS";
    echo
else
    echo "  - No InfiniBand"
fi

dnl if test -z "${HAVE_GENGETOPT_TRUE}"; then
dnl         echo "  - GENGETOPT";
dnl         echo "      - GENGETOPT_CFLAGS = $GENGETOPT_CFLAGS";
dnl         echo "      - GENGETOPT_CPPFLAGS = $GENGETOPT_CPPFLAGS";
dnl         echo "      - GENGETOPT_LDFLAGS = $GENGETOPT_LDFLAGS";
dnl         echo "      - GENGETOPT_LIBS = $GENGETOPT_LIBS";
dnl         echo
dnl else
dnl         echo "  - No Gengetopt"
dnl fi

if test "x${datatap}" != "xdisable"; then
    echo "  - DataTap with interconnect: $datatap";
    echo "      - DT_CFLAGS     = $DT_CFLAGS";
    echo "      - DT_CPPFLAGS   = $DT_CPPFLAGS";
    echo "      - DT_LDFLAGS    = $DT_LDFLAGS";
    echo "      - DT_LIBS       = $DT_LIBS";
    echo
else
    echo "  - No DataTap transport";
fi

if test -z "${HAVE_DART_TRUE}"; then
    echo "  - DART transport method";
    echo "      - DART_CFLAGS = $DART_CFLAGS";
    echo "      - DART_CPPFLAGS = $DART_CPPFLAGS";
    echo "      - DART_LDFLAGS = $DART_LDFLAGS";
    echo "      - DART_LIBS = $DART_LIBS";
    echo
else
    echo "  - No DART to build DART transport method"
fi

if test -z "${HAVE_DIMES_TRUE}"; then
    echo "  - DIMES transport method";
    echo "      - DIMES_CFLAGS = $DIMES_CFLAGS";
    echo "      - DIMES_CPPFLAGS = $DIMES_CPPFLAGS";
    echo "      - DIMES_LDFLAGS = $DIMES_LDFLAGS";
    echo "      - DIMES_LIBS = $DIMES_LIBS";
    echo
else
    echo "  - No DIMES to build DIMES transport method"
fi

if test -z "${HAVE_PHDF5_TRUE}"; then
    echo "  - PHDF5";
    echo "      - PHDF5_CFLAGS = $PHDF5_CFLAGS";
    echo "      - PHDF5_CPPFLAGS = $PHDF5_CPPFLAGS";
    echo "      - PHDF5_LDFLAGS = $PHDF5_LDFLAGS";
    echo "      - PHDF5_LIBS = $PHDF5_LIBS";
    echo
else
    echo "  - No PHDF5"
fi

if test -z "${HAVE_HDF5_TRUE}"; then
    echo "  - HDF5";
    echo "      - HDF5_CFLAGS = $HDF5_CFLAGS";
    echo "      - HDF5_CPPFLAGS = $HDF5_CPPFLAGS";
    echo "      - HDF5_LDFLAGS = $HDF5_LDFLAGS";
    echo "      - HDF5_LIBS = $HDF5_LIBS";
    echo
else
    echo "  - No HDF5"
fi

if test -z "${HAVE_NC4PAR_TRUE}"; then
    echo "  - NetCDF 4 Parallel";
    echo "      - NC4PAR_CFLAGS = $NC4PAR_CFLAGS";
    echo "      - NC4PAR_CPPFLAGS = $NC4PAR_CPPFLAGS";
    echo "      - NC4PAR_LDFLAGS = $NC4PAR_LDFLAGS";
    echo "      - NC4PAR_LIBS = $NC4PAR_LIBS";
    echo
else
    echo "  - No NETCDF 4 Parallel"
fi

if test -z "${HAVE_NETCDF_TRUE}"; then
    echo "  - NetCDF";
    echo "      - NETCDF_CFLAGS = $NETCDF_CFLAGS";
    echo "      - NETCDF_CPPFLAGS = $NETCDF_CPPFLAGS";
    echo "      - NETCDF_LDFLAGS = $NETCDF_LDFLAGS";
    echo "      - NETCDF_LIBS = $NETCDF_LIBS";
    echo
else
    echo "  - No NETCDF"
fi

if test -z "${HAVE_NSSI_TRUE}"; then
    echo "  - NSSI";
    echo "      - NSSI_CFLAGS = $NSSI_CFLAGS";
    echo "      - NSSI_CPPFLAGS = $NSSI_CPPFLAGS";
    echo "      - NSSI_LDFLAGS = $NSSI_LDFLAGS";
    echo "      - NSSI_LIBS = $NSSI_LIBS";
    echo
else
    echo "  - No NSSI"
fi

echo "  - SVN version = ${SVNVERSION}"
echo "  - Install path = ${prefix}"

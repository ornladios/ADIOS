#
#
# AC_PROG_DATASPACES
#
# Test for DATASPACES installation
# and set $DATASPACES to the correct value.
#
#
dnl @synopsis AC_PROG_DATASPACES
dnl
dnl This macro test if DATASPACES is installed. If DATASPACES
dnl is installed, it set $DATASPACES to the right value
dnl
dnl @version 1.0
dnl @author Norbert Podhorszki, pnorbert@ornl.gov
dnl
AC_DEFUN([AC_PROG_DATASPACES],[

AC_REQUIRE([AC_INFINIBAND])
dnl AC_REQUIRE([AC_PORTALS])
dnl AC_REQUIRE([AC_GNI])


AM_CONDITIONAL(HAVE_DATASPACES,true)

AC_ARG_WITH(dataspaces,
        [AS_HELP_STRING([--with-dataspaces=DIR],
           [Build the DATASPACES transport method. Point to the DATASPACES installation.])],
        [DATASPACES_LIBS="-L$withval/lib -ldspaces -ldscommon -ldart";
         DATASPACES_CPPFLAGS="-I$withval/include";
         DSPACES_CONFIG="$withval/bin/dspaces_config";
        ],
        [with_dataspaces=check])

dnl If --without-dataspaces was given set HAVE_DATASPACES to false and do nothing more
if test "x$with_dataspaces" == "xno"; then

   AM_CONDITIONAL(HAVE_DATASPACES,false)

elif test -z "${HAVE_MPI_FALSE}"; then

   AC_MSG_NOTICE([    DataSpaces does not build without MPI])
   AM_CONDITIONAL(HAVE_DATASPACES,false)

else

    dnl allow args --with-dataspaces incdir and --with-dataspaces-libdir
    AC_ARG_WITH(dataspaces-incdir,
                [  --with-dataspaces-incdir=<location of dataspaces includes>],
                [DATASPACES_INCDIR=$withval
                 with_dataspaces=detailed])
    
    AC_ARG_WITH(dataspaces-libdir,
                [  --with-dataspaces-libdir=<location of dataspaces library>],
                [DATASPACES_LIBDIR=$withval
                 with_dataspaces=detailed])
    
    
    dnl If we know DATASPACES_DIR, then we can know DATASPACES_INCDIR.
    dnl We don't overwrite DATASPACES_INCDIR.
    if test -n "${DATASPACES_DIR}" -a -z "${DATASPACES_INCDIR}"; then
            DATASPACES_INCDIR="${DATASPACES_DIR}/include";
    dnl We may have DATASPACES denoting the dir (e.g. on ewok BUT on franklin it contains all flags)
    elif test -n "${DATASPACES}" -a -d "${DATASPACES}"; then
            DATASPACES_INCDIR="${DATASPACES}/include"
    fi

    dnl If we know DATASPACES_DIR, then we can know DATASPACES_LIBDIR.
    dnl We don't overwrite DATASPACES_LIBDIR.
    if test -n "${DATASPACES_DIR}" -a -z "${DATASPACES_LIBDIR}"; then
            DATASPACES_LIBDIR="${DATASPACES_DIR}/lib";
    dnl We may have DATASPACES denoting the dir (e.g. on ewok BUT on franklin it contains all flags)
    elif test -n "${DATASPACES}" -a -d "${DATASPACES}"; then
            DATASPACES_LIBDIR="${DATASPACES}/lib"
    fi

    dnl 1.6.2 has a config file
    if test -x "${DSPACES_CONFIG}"; then
        DATASPACES_LIBS=$(${DSPACES_CONFIG} -l)
        DATASPACES_CPPFLAGS=$(${DSPACES_CONFIG} -c)
    else

        dnl Add "-I" to DATASPACES_INCDIR.
        if test -n "${DATASPACES_INCDIR}"; then
                DATASPACES_CPPFLAGS="-I${DATASPACES_INCDIR}"
        fi

        dnl Add "-L" to DATASPACES_LIBDIR.
        if test -n "${DATASPACES_LIBDIR}"; then
                DATASPACES_LIBS="-L${DATASPACES_LIBDIR} -ldspaces -ldscommon -ldart"
        fi
    fi

    save_CC="$CC"
    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    dnl LDFLAGS="$LDFLAGS $DATASPACES_LDFLAGS"
    LIBS="$DATASPACES_LIBS $LIBS"
    CPPFLAGS="$CPPFLAGS $DATASPACES_CPPFLAGS"
    CC="$MPICC"
    
    if test -z "${HAVE_DATASPACES_TRUE}"; then
            AC_CHECK_HEADERS(dataspaces.h,
                    ,
                    [AM_CONDITIONAL(HAVE_DATASPACES,false)])
    fi
    
    if test -z "${HAVE_DATASPACES_TRUE}"; then
        # Check for the DataSpaces library and headers
        AC_MSG_CHECKING([if dataspaces code can be compiled and linked])
        dnl if test "x${ac_portals_lib_ok}" == "xyes"; then 
        dnl elif test "x${ac_infiniband_lib_ok}" == "xyes"; then 
        dnl elif test -z "${HAVE_CRAY_PMI_TRUE}" -a -z "${HAVE_CRAY_UGNI_TRUE}"; then 
        dnl elif test "x${ac_dcmf_lib_ok}" == "xyes"; then
        dnl elif test "x${ac_pami_lib_ok}" == "xyes"; then
        dnl else
        dnl fi
            AC_TRY_LINK([#include "dataspaces.h"],
                    [int err; err = dspaces_init(1,1,0,"");],
                    [AC_MSG_RESULT(yes)],
                    [AC_MSG_RESULT(no)
                     AM_CONDITIONAL(HAVE_DATASPACES,false)])
    fi

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    CC="$save_CC"
    
    AC_SUBST(DATASPACES_LIBS)
    AC_SUBST(DATASPACES_LDFLAGS)
    AC_SUBST(DATASPACES_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_DATASPACES_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_DATASPACES,1,[Define if you have the DATASPACES.])],[$1])
            :
    else
            $2
            :
    fi

fi
])dnl AC_DATASPACES

#
#
# AC_PROG_DIMES
#
# Test for DIMES installation
# and set $DIMES to the correct value.
#
#

dnl @synopsis AC_PROG_DIMES
dnl
dnl This macro test if DIMES is installed. If DIMES
dnl is installed, it set $DIMES to the right value
dnl
dnl @version 1.0
dnl @author Fan Zhang,  zhangfan@cac.rutgers.edu
dnl @author Norbert Podhorszki, pnorbert@ornl.gov
dnl
AC_DEFUN([AC_PROG_DIMES],[

AM_CONDITIONAL(HAVE_DIMES,true)

AC_ARG_WITH(dimes,
        [AS_HELP_STRING([--with-dimes=DIR],
           [Build the DIMES transport method. Point to the DIMES installation.])],
        [DIMES_LDFLAGS="-L$withval/lib";
         DIMES_CPPFLAGS="-I$withval/include";],
        [with_dimes=check])

dnl If --without-dimes was given set HAVE_DIMES to false and do nothing more
if test "x$with_dimes" == "xno"; then

   AM_CONDITIONAL(HAVE_DIMES,false)

elif test -z "${HAVE_MPI_FALSE}"; then

   AC_MSG_NOTICE([    Dimes does not build without MPI])
   AM_CONDITIONAL(HAVE_DIMES,false)

else

    dnl allow args --with-dimes incdir and --with-dimes-libdir
    AC_ARG_WITH(dimes-incdir,
                [  --with-dimes-incdir=<location of dimes includes>],
                [DIMES_INCDIR=$withval
                 with_dimes=detailed])
    
    AC_ARG_WITH(dimes-libdir,
                [  --with-dimes-libdir=<location of dimes library>],
                [DIMES_LIBDIR=$withval
                 with_dimes=detailed])
    
    
    dnl If we know DIMES_DIR, then we can know DIMES_INCDIR.
    dnl We don't overwrite DIMES_INCDIR.
    if test -n "${DIMES_DIR}" -a -z "${DIMES_INCDIR}"; then
            DIMES_INCDIR="${DIMES_DIR}/include";
    dnl We may have DIMES denoting the dir (e.g. on ewok BUT on franklin it contains all flags)
    elif test -n "${DIMES}" -a -d "${DIMES}"; then
            DIMES_INCDIR="${DIMES}/include"
    fi

    dnl If we know DIMES_DIR, then we can know DIMES_LIBDIR.
    dnl We don't overwrite DIMES_LIBDIR.
    if test -n "${DIMES_DIR}" -a -z "${DIMES_LIBDIR}"; then
            DIMES_LIBDIR="${DIMES_DIR}/lib";
    dnl We may have DIMES denoting the dir (e.g. on ewok BUT on franklin it contains all flags)
    elif test -n "${DIMES}" -a -d "${DIMES}"; then
            DIMES_LIBDIR="${DIMES}/lib"
    fi

    dnl Add "-I" to DIMES_INCDIR.
    if test -n "${DIMES_INCDIR}"; then
            DIMES_CPPFLAGS="-I${DIMES_INCDIR}"
    else
            ac_dimes_ok=no
    fi

    dnl Add "-L" to DIMES_LIBDIR.
    if test -n "${DIMES_LIBDIR}"; then
            DIMES_LDFLAGS="-L${DIMES_LIBDIR}"
    else
            ac_dimes_ok=no
    fi

    save_CC="$CC"
    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    if test "x${ac_portals_lib_ok}" == "xyes"; then 
        dnl LIBS="$LIBS -ldart2 -lspaces"
	echo "DIMES currently NOT supported for Cray Portals!"
	AM_CONDITIONAL(HAVE_DIMES, false)
    else
        LIBS="-ldspaces -ldscommon -ldart $LIBS"
    fi
    LDFLAGS="$LDFLAGS $DIMES_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $DIMES_CPPFLAGS"
    CC="$MPICC"
    
    if test -z "${HAVE_DIMES_TRUE}"; then
            AC_CHECK_HEADERS(dimes_interface.h,
                    ,
                    [AM_CONDITIONAL(HAVE_DIMES,false)])
    fi
    
    if test -z "${HAVE_DIMES_TRUE}"; then
        AC_MSG_CHECKING([if dimes code can be compiled and linked])
        # Check for the DataSpaces/DIMES library and headers
        dnl if test -z "${HAVE_CRAY_PMI_TRUE}" -a -z "${HAVE_CRAY_UGNI_TRUE}"; then 
	    dnl elif test "x${ac_infiniband_lib_ok}" == "xyes"; then 
	    dnl elif test "x${ac_dcmf_lib_ok}" == "xyes"; then
	    dnl elif test "x${ac_pami_lib_ok}" == "xyes"; then
           AC_TRY_LINK([#include "dimes_interface.h"],
                    [int err; dimes_put_sync_all();],
                    [AC_MSG_RESULT(yes)
                     DIMES_LIBS="-ldspaces -ldscommon -ldart"],
                    [AC_MSG_RESULT(no)
                     AM_CONDITIONAL(HAVE_DIMES,false)])
        dnl else
            dnl AM_CONDITIONAL(HAVE_DIMES,false)
        dnl fi
    fi
 
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    CC="$save_CC"
    
    AC_SUBST(DIMES_LIBS)
    AC_SUBST(DIMES_LDFLAGS)
    AC_SUBST(DIMES_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_DIMES_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_DIMES,1,[Define if you have the DIMES.])],[$1])
            :
    else
            $2
            :
    fi

fi
])dnl AC_DIMES

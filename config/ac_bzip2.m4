#
#
# AC_BZIP2
#
#
#
dnl @synopsis AC_BZIP2
dnl
dnl This macro test if BZIP2 is to be used.
dnl Use in C code:
dnl     #ifdef BZIP2
dnl     #include "bzlib.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Zhenhuan (Steve) Gong
dnl
AC_DEFUN([AC_BZIP2],[

AC_MSG_NOTICE([=== checking for BZIP2 ===])

AM_CONDITIONAL(HAVE_BZIP2,true)

AC_ARG_WITH(bzip2,
        [  --with-bzip2=DIR      Location of BZIP2 library],
        [BZIP2_LDFLAGS="-L$withval/lib";
         BZIP2_LIBS="-lbz2";
         BZIP2_CPPFLAGS="-I$withval/include";],
        [with_bzip2=no])

if test "x$with_bzip2" == "xno"; then

   AM_CONDITIONAL(HAVE_BZIP2,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lbz2"
    LDFLAGS="$LDFLAGS $BZIP2_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $BZIP2_CPPFLAGS"

    if test -z "${HAVE_BZIP2_TRUE}"; then
           AC_CHECK_HEADERS(bzlib.h,
                   ,
                   [AM_CONDITIONAL(HAVE_BZIP2,false)])
    fi

    # Check for the BZIP2 library and headers
    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [BZIP2_LIBS="-lbz2"],
    dnl        [AM_CONDITIONAL(HAVE_BZIP2,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(BZIP2_LIBS)
    AC_SUBST(BZIP2_LDFLAGS)
    AC_SUBST(BZIP2_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_BZIP2_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_BZIP2,1,[Define if you have BZIP2.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_BZIP2

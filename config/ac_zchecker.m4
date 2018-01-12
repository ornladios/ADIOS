#
#
# AC_ZCHECKER
#
#
#
dnl @synopsis AC_ZCHECKER
dnl
dnl This macro test if Z-Checker is to be used.
dnl Use in C code:
dnl     #ifdef HAVE_ZCHECKER
dnl     #include "zc.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Zhenhuan (Steve) Gong
dnl
AC_DEFUN([AC_ZCHECKER],[

AC_MSG_NOTICE([=== checking for Z-Checker ===])

AM_CONDITIONAL(HAVE_ZCHECKER,true)

AC_ARG_WITH(zchecker,
        [  --with-zchecker=DIR      Location of Z-Checker library],
        [ZCHECKER_LDFLAGS="-L$withval/lib";
         ZCHECKER_LIBS="-lzc";
         ZCHECKER_CPPFLAGS="-I$withval/include";],
        [with_ZCHECKER=no])

if test "x$with_ZCHECKER" == "xno"; then

   AM_CONDITIONAL(HAVE_ZCHECKER,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lzc"
    LDFLAGS="$LDFLAGS $ZCHECKER_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $ZCHECKER_CPPFLAGS"

    dnl if test -z "${HAVE_ZCHECKER_TRUE}"; then
    dnl        AC_CHECK_HEADERS(zc.h,
    dnl                ,
    dnl                [AM_CONDITIONAL(HAVE_ZCHECKER,false)])
    dnl fi

    # Check for the Z-Checker library and headers
    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [ZCHECKER_LIBS="-lzc"],
    dnl        [AM_CONDITIONAL(HAVE_ZCHECKER,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(ZCHECKER_LIBS)
    AC_SUBST(ZCHECKER_LDFLAGS)
    AC_SUBST(ZCHECKER_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_ZCHECKER_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_ZCHECKER,1,[Define if you have Z-Checker.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_ZCHECKER

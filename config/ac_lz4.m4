#
#
# AC_LZ4
#
#
#
dnl @synopsis AC_LZ4
dnl
dnl This macro test if LZ4 is to be used.
dnl Use in C code:
dnl     #ifdef LZ4
dnl     #include "lz4.h"
dnl     #endif
dnl
dnl @version 2.0
dnl @author Zhenhuan (Steve) Gong
dnl @author Norbert Podhorszki
dnl @author Rene Widera
dnl
AC_DEFUN([AC_LZ4],[

AC_MSG_NOTICE([=== checking for LZ4 ===])

AM_CONDITIONAL(HAVE_LZ4,true)

AC_ARG_WITH(lz4,
        [  --with-lz4=DIR      Location of LZ4 library],
        [:], [with_lz4=no])

if test "x$with_lz4" == "xno"; then

   AM_CONDITIONAL(HAVE_LZ4,false)

else

    save_CFLAGS="$CFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"

    if test "x$with_lz4" == "xyes"; then
        dnl No path given
        LZ4_CFLAGS=""
        LZ4_CPPFLAGS=""
        LZ4_LDFLAGS=""
        LZ4_LIBS="-llz4"
    else
        dnl Path given, first try path/lib
        LZ4_CFLAGS="-I$withval/include"
        LZ4_CPPFLAGS="-I$withval/include"
        LZ4_LDFLAGS="-L$withval/lib"
        LZ4_LIBS="-llz4"
    fi

    LIBS="$LIBS $LZ4_LIBS"
    LDFLAGS="$LDFLAGS $LZ4_LDFLAGS"
    CFLAGS="$CFLAGS $LZ4_CFLAGS"
    CPPFLAGS="$CPPFLAGS $LZ4_CPPFLAGS"

    dnl Find header file first
    AC_CHECK_HEADERS(lz4.h,
              ,
              [AM_CONDITIONAL(HAVE_LZ4,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CFLAGS="$save_CFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(LZ4_LIBS)
    AC_SUBST(LZ4_LDFLAGS)
    AC_SUBST(LZ4_CFLAGS)
    AC_SUBST(LZ4_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_LZ4_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_LZ4,1,[Define if you have LZ4.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_LZ4

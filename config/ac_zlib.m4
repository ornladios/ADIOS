#
#
# AC_ZLIB
#
#
#
dnl @synopsis AC_ZLIB
dnl
dnl This macro test if ZLIB is to be used.
dnl Use in C code:
dnl     #ifdef ZLIB
dnl     #include "zlib.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Zhenhuan (Steve) Gong
dnl
AC_DEFUN([AC_ZLIB],[

AC_MSG_NOTICE([=== checking for ZLIB ===])

AM_CONDITIONAL(HAVE_ZLIB,true)

AC_ARG_WITH(zlib,
        [  --with-zlib=DIR      Location of ZLIB library],
        [ZLIB_LDFLAGS="-L$withval/lib -L$withval/lib64";
         ZLIB_LIBS="-lz";
         ZLIB_CPPFLAGS="-I$withval/include";],
        [with_zlib=no])

if test "x$with_zlib" == "xno"; then

   AM_CONDITIONAL(HAVE_ZLIB,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lz"
    LDFLAGS="$LDFLAGS $ZLIB_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $ZLIB_CPPFLAGS"

    if test -z "${HAVE_ZLIB_TRUE}"; then
           AC_CHECK_HEADERS(zlib.h,
                   ,
                   [AM_CONDITIONAL(HAVE_ZLIB,false)])
    fi

    # Check for the ZLIB library and headers
    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [ZLIB_LIBS="-lz"],
    dnl        [AM_CONDITIONAL(HAVE_ZLIB,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(ZLIB_LIBS)
    AC_SUBST(ZLIB_LDFLAGS)
    AC_SUBST(ZLIB_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_ZLIB_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_ZLIB,1,[Define if you have ZLIB.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_ZLIB

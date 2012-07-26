#
#
# AC_COMPRESS
#
#
#
dnl @synopsis AC_COMPRESS
dnl
dnl This macro test if COMPRESS is to be used.
dnl Use in C code:
dnl     #ifdef COMPRESS
dnl     #include "compress.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Zhenhuan (Steve) Gong
dnl
AC_DEFUN([AC_COMPRESS],[

AC_MSG_NOTICE([=== checking for COMPRESS ===])

AM_CONDITIONAL(HAVE_COMPRESS,true)

AC_ARG_WITH(compress,
        [  --with-compress=DIR      Location of COMPRESS library],
        [COMPRESS_LDFLAGS="-L$withval/lib";
         COMPRESS_LIBS="-lcompress -lbspline -lbz2 -lz";
         COMPRESS_CPPFLAGS="-I$withval/include";],
        [with_compress=no])

if test "x$with_compress" == "xno"; then

   AM_CONDITIONAL(HAVE_COMPRESS,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lcompress -lbspline -lbz2 -lz"
    LDFLAGS="$LDFLAGS $COMPRESS_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $COMPRESS_CPPFLAGS"

    dnl if test -z "${HAVE_COMPRESS_TRUE}"; then
    dnl        AC_CHECK_HEADERS(compress.h,
    dnl                ,
    dnl                [AM_CONDITIONAL(HAVE_COMPRESS,false)])
    dnl fi

    # Check for the COMPRESS library and headers
    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [COMPRESS_LIBS="-lcompress -lbspline -lbz2 -lz"],
    dnl        [AM_CONDITIONAL(HAVE_COMPRESS,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(COMPRESS_LIBS)
    AC_SUBST(COMPRESS_LDFLAGS)
    AC_SUBST(COMPRESS_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_COMPRESS_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_COMPRESS,1,[Define if you have COMPRESS.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_COMPRESS

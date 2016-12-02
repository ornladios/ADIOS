#
#
# AC_SZ
#
#
#
dnl @synopsis AC_SZ
dnl
dnl This macro test if SZ is to be used.
dnl Use in C code:
dnl     #ifdef SZ
dnl     #include "sz.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Zhenhuan (Steve) Gong
dnl
AC_DEFUN([AC_SZ],[

AC_MSG_NOTICE([=== checking for SZ ===])

AM_CONDITIONAL(HAVE_SZ,true)

AC_ARG_WITH(sz,
        [  --with-sz=DIR      Location of SZ library],
        [SZ_LDFLAGS="-L$withval/lib";
         SZ_LIBS="-lsz -lzlib";
         SZ_CPPFLAGS="-I$withval/include";],
        [with_sz=no])

if test "x$with_sz" == "xno"; then

   AM_CONDITIONAL(HAVE_SZ,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lsz -lzlib"
    LDFLAGS="$LDFLAGS $SZ_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $SZ_CPPFLAGS"

    if test -z "${HAVE_SZ_TRUE}"; then
           AC_CHECK_HEADERS(sz.h,
                   ,
                   [AM_CONDITIONAL(HAVE_SZ,false)])
    fi

    # Check for the SZ library and headers
    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [SZ_LIBS="-lsz"],
    dnl        [AM_CONDITIONAL(HAVE_SZ,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(SZ_LIBS)
    AC_SUBST(SZ_LDFLAGS)
    AC_SUBST(SZ_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_SZ_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_SZ,1,[Define if you have SZ.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_SZ

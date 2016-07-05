#
#
# AC_ZFP
#
#
#
dnl @synopsis AC_ZFP
dnl
dnl This macro tests if ZFP is to be used.
dnl Use in C code:
dnl     #ifdef ZFP
dnl     #include "zfp.h"
dnl     #endif
dnl
dnl @version 0.1
dnl @author Eric Suchyta
dnl

AC_DEFUN([AC_ZFP],[

AC_MSG_NOTICE([=== checking for ZFP ===])
AM_CONDITIONAL(HAVE_ZFP,true)

AC_ARG_WITH(zfp,
        [  --with-zfp=DIR      Location of ZFP library],
        [ZFP_LDFLAGS="-L$withval/lib";
         ZFP_LIBS="-lzfp";
         ZFP_CPPFLAGS="-I$withval/include";],
        [with_zfp=no])

if test "x$with_zfp" == "xno"; then
   AM_CONDITIONAL(HAVE_ZFP,false)

else
    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lzfp"
    LDFLAGS="$LDFLAGS $ZFP_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $ZFP_CPPFLAGS"

    if test -z "${HAVE_ZFP_TRUE}"; then
           AC_CHECK_HEADERS(zfp.h, , [AM_CONDITIONAL(HAVE_ZFP,false)])
    fi

    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [SZIP_LIBS="-lsz"],
    dnl        [AM_CONDITIONAL(HAVE_SZIP,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(ZFP_LIBS)
    AC_SUBST(ZFP_LDFLAGS)
    AC_SUBST(ZFP_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_ZFP_TRUE}"; then
            ifelse([$1], , [AC_DEFINE(HAVE_ZFP,1,[Define if you have ZFP.])], [$1])
            :
    else
            $2
            :
    fi
fi
])
dnl AC_ZFP

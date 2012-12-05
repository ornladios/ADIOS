#
#
# AC_MLOC
#
#
#
dnl @synopsis AC_MLOC
dnl
dnl This macro test if MLOC is to be used.
dnl Use in C code:
dnl     #ifdef MLOC
dnl     #include "mloc.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Zhenhuan (Steve) Gong
dnl
AC_DEFUN([AC_MLOC],[

AC_MSG_NOTICE([=== checking for MLOC ===])

AM_CONDITIONAL(HAVE_MLOC,true)

AC_ARG_WITH(mloc,
        [  --with-mloc=DIR      Location of MLOC library],
        [MLOC_LDFLAGS="-L$withval/lib";
         MLOC_LIBS="-lmloc";
         MLOC_CPPFLAGS="-I$withval/include";],
        [with_mloc=no])

if test "x$with_mloc" == "xno"; then

   AM_CONDITIONAL(HAVE_MLOC,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lmloc"
    LDFLAGS="$LDFLAGS $MLOC_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $MLOC_CPPFLAGS"

    dnl if test -z "${HAVE_MLOC_TRUE}"; then
    dnl        AC_CHECK_HEADERS(mloc.h,
    dnl                ,
    dnl                [AM_CONDITIONAL(HAVE_MLOC,false)])
    dnl fi

    # Check for the MLOC library and headers
    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [MLOC_LIBS="-lmloc"],
    dnl        [AM_CONDITIONAL(HAVE_MLOC,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(MLOC_LIBS)
    AC_SUBST(MLOC_LDFLAGS)
    AC_SUBST(MLOC_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_MLOC_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_MLOC,1,[Define if you have MLOC.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_MLOC

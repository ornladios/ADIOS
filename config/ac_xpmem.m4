#
#
# AC_XPMEM
#
#
#
dnl @synopsis AC_XPMEM
dnl
dnl This macro test if XPMEM is to be used.
dnl Use in C code:
dnl     #ifdef XPMEM
dnl     #include "szlib.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Zhenhuan (Steve) Gong
dnl
AC_DEFUN([AC_XPMEM],[

AC_MSG_NOTICE([=== checking for XPMEM ===])

AM_CONDITIONAL(HAVE_XPMEM,true)

AC_ARG_WITH(xpmem,
        [  --with-xpmem=DIR      Location of XPMEM library],
        [XPMEM_LDFLAGS="-L$withval/lib";
         XPMEM_LIBS="-lxpmem";
         XPMEM_CPPFLAGS="-I$withval/include";],
        [with_xpmem=no])

if test "x$with_xpmem" == "xno"; then

   AM_CONDITIONAL(HAVE_XPMEM,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lxpmem"
    LDFLAGS="$LDFLAGS $XPMEM_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $XPMEM_CPPFLAGS"

    if test -z "${HAVE_XPMEM_TRUE}"; then
           AC_CHECK_HEADERS(xpmem.h,
                   ,
                   [AM_CONDITIONAL(HAVE_XPMEM,false)])
    fi

    # Check for the XPMEM library and headers
    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [XPMEM_LIBS="-l"],
    dnl        [AM_CONDITIONAL(HAVE_XPMEM,false)])

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(XPMEM_LIBS)
    AC_SUBST(XPMEM_LDFLAGS)
    AC_SUBST(XPMEM_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_XPMEM_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_XPMEM,1,[Define if you have XPMEM.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_XPMEM

#
#
# AC_MGARD
#
#
#
dnl @synopsis AC_MGARD
dnl
dnl This macro tests if MGARD is to be used.
dnl Use in C code:
dnl     #ifdef MGARD
dnl     #include "mgard_capi.h"
dnl     #endif
dnl
dnl @version 0.1
dnl @author Jong Choi
dnl

AC_DEFUN([AC_MGARD],[

AC_MSG_NOTICE([=== checking for MGARD ===])

AM_CONDITIONAL(HAVE_MGARD,true)

AC_ARG_WITH(mgard,
        [  --with-mgard=DIR      Location of MGARD library],
        [],
        [with_mgard=no])


if test "x$with_mgard" == "xno"; then

    AM_CONDITIONAL(HAVE_MGARD,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"

    if test -z "${MGARD_LIBS}"; then
            MGARD_LIBS="-lmgard -lstdc++ -lz -lm"
    fi
    
    if test -z "${MGARD_LDFLAGS}"; then
            MGARD_LDFLAGS="-L$with_mgard/lib"
    fi
    
    if test -z "${MGARD_CPPFLAGS}"; then
            MGARD_CPPFLAGS="-I$with_mgard/include"
    fi

    LIBS="$LIBS $MGARD_LIBS"
    LDFLAGS="$LDFLAGS $MGARD_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $MGARD_CPPFLAGS"

    if test -z "${HAVE_MGARD_TRUE}"; then
           AC_CHECK_HEADERS(mgard_capi.h,,[AM_CONDITIONAL(HAVE_MGARD,false)])
    fi

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(MGARD_LIBS)
    AC_SUBST(MGARD_LDFLAGS)
    AC_SUBST(MGARD_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_MGARD_TRUE}"; then
            ifelse([$1], , [AC_DEFINE(HAVE_MGARD,1,[Define if you have MGARD.])], [$1])
            :
    else
            $2
            :
    fi
fi
])
dnl AC_MGARD

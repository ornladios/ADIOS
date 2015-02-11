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
dnl     #include <xpmem.h>
dnl     #endif
dnl
dnl @version 2.0
dnl @author Zhenhuan (Steve) Gong
dnl dnl @author Norbert Podhorszki
dnl dnl @modified Hasan Abbasi to work with XPMEM

AC_DEFUN([AC_XPMEM],[

AC_MSG_NOTICE([=== checking for XPMEM ===])

AM_CONDITIONAL(HAVE_XPMEM,true)

AC_ARG_WITH(xpmem,
        [  --with-xpmem=DIR      Location of XPMEM library],
        [:], [with_xpmem=no])

if test "x$with_xpmem" == "xno"; then

   AM_CONDITIONAL(HAVE_XPMEM,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"

    if test "x$with_xpmem" == "xyes"; then
        dnl No path given
        XPMEM_CPPFLAGS=""
        XPMEM_LDFLAGS=""
        XPMEM_LIBS="-lxpmem"
    else
        dnl Path given, first try path/lib
        XPMEM_CPPFLAGS="-I$withval/include"
        XPMEM_LDFLAGS="-L$withval/lib"
        XPMEM_LIBS="-lxpmem"
    fi

    LIBS="$LIBS $XPMEM_LIBS"
    LDFLAGS="$LDFLAGS $XPMEM_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $XPMEM_CPPFLAGS"

    if test -z "${HAVE_XPMEM_TRUE}"; then
           AC_CHECK_HEADERS(xpmem.h,
                   ,
                   [AM_CONDITIONAL(HAVE_XPMEM,false)])
    fi

    if test -z "${HAVE_XPMEM_TRUE}"; then
        dnl Try to link an example now
        AC_MSG_CHECKING([if xpmem code can be linked with $XPMEM_LDFLAGS])
        AC_TRY_LINK(
            [#include <stdlib.h>
             #include <xpmem.h>],
            [return xpmem_init();],
            [AC_MSG_RESULT(yes)],
            [AM_CONDITIONAL(HAVE_XPMEM,false)
             AC_MSG_RESULT(no)
            ])

        dnl If linking above failed, one reason might be that we looked in lib/
        dnl instead of lib64/
        if test -z "${HAVE_XPMEM_FALSE}"; then
            XPMEM_LDFLAGS="-L$withval/lib64"
            LDFLAGS="$LDFLAGS $XPMEM_LDFLAGS"
            AC_MSG_CHECKING([if xpmem code can be linked with $XPMEM_LDFLAGS])
            AC_TRY_LINK(
            [#include <stdlib.h>
             #include <xpmem.h>],
            [return xpmem_init();],
            [AC_MSG_RESULT(yes)],
            [AM_CONDITIONAL(HAVE_XPMEM,false)
             AC_MSG_RESULT(no)
            ])
        fi
    fi
 

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

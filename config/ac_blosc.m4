#
#
# AC_BLOSC
#
#
#
dnl @synopsis AC_BLOSC
dnl
dnl This macro test if blosc (for Titan only) is to be used. 
dnl Use in C code:
dnl     #ifdef BLOSC
dnl     #include "blosc.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Qing Liu, ORNL
dnl
AC_DEFUN([AC_BLOSC],[

AC_MSG_NOTICE([=== checking for BLOSC ===])

AM_CONDITIONAL(HAVE_BLOSC,true)

AC_ARG_WITH(blosc,
        [  --with-blosc=DIR      Location of BLOSC library],
        [BLOSC_LDFLAGS="-L$withval";
         BLOSC_LIBS="-L/ccs/proj/e2e/qliu/blosc/lib -lblosc -pthread /opt/gcc/6.3.0/snos/lib/../lib64/libstdc++.a";
         BLOSC_CPPFLAGS="-I$withval/include";],
        [with_blosc=no])

if test "x$with_blosc" == "xno"; then

   AM_CONDITIONAL(HAVE_BLOSC,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lblosc"
    LDFLAGS="$LDFLAGS $FDR_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $FDR_CPPFLAGS"
    
    # Check for the BLOSC library and headers
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    
    AC_SUBST(BLOSC_LIBS)
    AC_SUBST(BLOSC_LDFLAGS)
    AC_SUBST(BLOSC_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_BLOSC_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_BLOSC,1,[Define if you have BLOSC.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_BLOSC

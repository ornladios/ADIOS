#
#
# AC_MGARD
#
#
#
dnl @synopsis AC_MGARD
dnl
dnl This macro test if mgard (for Titan only) is to be used. 
dnl Use in C code:
dnl     #ifdef MGARD
dnl     #include "mgard.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Qing Liu, ORNL
dnl
AC_DEFUN([AC_MGARD],[

AC_MSG_NOTICE([=== checking for MGARD ===])

AM_CONDITIONAL(HAVE_MGARD,true)

AC_ARG_WITH(mgard,
        [  --with-mgard=DIR      Location of MGARD library],
        [MGARD_LDFLAGS="-L$withval";
         MGARD_LIBS="-lmgard -lz -lm /opt/gcc/6.3.0/snos/lib/../lib64/libstdc++.a -L/ccs/proj/e2e/qliu/blosc/lib -lblosc -pthread /opt/gcc/6.3.0/snos/lib/../lib64/libstdc++.a";
         MGARD_CPPFLAGS="-I$withval/include";],
        [with_mgard=no])

if test "x$with_mgard" == "xno"; then

   AM_CONDITIONAL(HAVE_MGARD,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lmgard"
    LDFLAGS="$LDFLAGS $FDR_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $FDR_CPPFLAGS"
    
    # Check for the MGARD library and headers
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    
    AC_SUBST(MGARD_LIBS)
    AC_SUBST(MGARD_LDFLAGS)
    AC_SUBST(MGARD_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_MGARD_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_MGARD,1,[Define if you have MGARD.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_MGARD
